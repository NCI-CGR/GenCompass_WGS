#!/usr/bin/env python3

# usage: cromreporter.py [-h] [-n NAME] [-w [WORKFLOW_LOGS ...]] [-wl WORKFLOW_LOG_LIST] [-wr WORKFLOW_ROOT] [-m {slurm,local}] [-s | --save-individual-reports | --no-save-individual-reports] [-a [{runtime,memory,gantt} ...]] [-odir ODIR]

# Generates a report of all subjobs

# options:
#   -h, --help            show this help message and exit
#   -n NAME, --name NAME
#   -w [WORKFLOW_LOGS ...], --workflow_logs [WORKFLOW_LOGS ...]
#                         workflow logs
#   -wl WORKFLOW_LOG_LIST, --workflow_log_list WORKFLOW_LOG_LIST
#                         List of workflow logs
#   -wr WORKFLOW_ROOT, --workflow_root WORKFLOW_ROOT
#                         Root directory of workflow logs. Directories up to 2 levels below are searched for logs
#   -m {slurm,local}, --mode {slurm,local}
#                         Job manager. If slurm then additional metrics are generated.
#   -s, --save-individual-reports, --no-save-individual-reports
#                         whether to save individual workflow reports (default: False)
#   -a [{runtime,memory,gantt} ...], --analysis [{runtime,memory,gantt} ...]
#                         Analysis to run in addition to the tsv report
#   -odir ODIR, --output_directory ODIR
#                         Output directory. Directory is created if it does not exist.

import pandas as pd
import re
import subprocess
from collections import defaultdict
from io import StringIO
import argparse
import os
from datetime import datetime
import logging
import sys
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
ch = logging.StreamHandler(sys.stderr)
ch.setLevel(logging.DEBUG)
logger.addHandler(ch)


class CromReporter:
    def __init__(self, logfile: str, workflow_name: str = 'workflow', mode: str = 'slurm'):
        """Generates report from Cromwell workflow logfile. 

        Args:
            logfile (str): filename of logfile
            workflow_name (str, optional): Name of the workflow. Defaults to 'workflow'.
            mode (str, optional): Job manager. Defaults to 'slurm'.
        """
        self.logfile: str = logfile
        self.workflow_name: str = workflow_name
     
        self.mode: str = mode
        if self.logfile.startswith('gs://'):
            bucket, blob_name = split_gcp_filename(self.logfile)
            self.blob = download_gcp_blob(bucket, blob_name)
        else:
            self.blob = None

        self.workflow_id: str = self.get_workflow_id()
        self.task_logs_dict : defaultdict = defaultdict(lambda: defaultdict(lambda: defaultdict(str)))
        self.report: pd.DataFrame = None

        self.read_logfile()

    
    def get_workflow_id(self):
        """Extracts the base workflow id of the logfile.
        Base workflow id is extracted from the first line of the logfile

        Returns:
            workflow_id: logfile base workflow_id
        """
        if self.blob is None:
            with open(self.logfile) as log:
                firstline = log.readline()
        else:
            with self.blob.open("r") as f:
                firstline = f.readline()
        logger.info(firstline)
        result = re.search('UUID\((.*)\)', firstline)
        if result is not None:
            workflow_id = result.group(1)
            logger.info(f"workflow id {workflow_id}")
            return workflow_id
        else:
            return None

    

    
    def read_logfile_line(self, line: str):
        """Reads the logfile line and updates the logfile dict if line contains 
        relevant information.

        Args:
            line (str): logfile line
        """
        if 'job id' in line:
            logger.info(f"job id in line \n {line}")
            result = re.search('UUID\((.*)\)+(.*)]: job id: (.*)', line)
            workflow = result.group(1)
            task_id = result.group(2)
            task_name = re.search(r'\.(.*)\:(.*):(.*)', task_id).group(1)
            job_id = result.group(3)
            logger.info(f"workflow: {workflow}\ttask id: {task_id}\tjob_id: {job_id}")
            self.task_logs_dict[workflow][task_id]['JobID'] = job_id
            self.task_logs_dict[workflow][task_id]['TaskName'] = task_name
        elif 'Status change' in line:
            result = re.search('UUID\((.*)\)+(.*)]: Status change from (.*) to (.*)', line)
            workflow = result.group(1)
            task = result.group(2)
            old_status = result.group(3)
            new_status = result.group(4)
            self.task_logs_dict[workflow][task]['Status'] = new_status
            if self.mode == 'local':
                if new_status == "Running":
                    split = line.split()
                    self.task_logs_dict[workflow][task]['Start'] = f'{split[0]} {split[1]}'
                if new_status in ("Done", "Failed", "Success"):
                    split = line.split()
                    self.task_logs_dict[workflow][task]['End'] = f'{split[0]} {split[1]}'

    
    def update_report_with_sacct_data(self):
        """Update the report with additional generated from slurm sacct command. 
        Additional columns capture CPU and memory performance.
        """
        
        job_ids = self.report['JobID'].unique()
        job_array = [['-j', job] for job in job_ids]
        job_array = [item for sublist in job_array for item in sublist]
        fields="JobID%30,JobName,User,Partition,NodeList,Elapsed,Timelimit,State,ExitCode,MaxRSS,AllocTRES%32,ReqMem,AllocCPUS,MinCPU,AveCPU,MaxDiskRead,MaxDiskWrite,Start,End"
        sacct_command = ["sacct", "-P", "-o", fields] + job_array
        
        sacct_output = subprocess.check_output(sacct_command)
        sacct_output = sacct_output.decode('utf-8')
        sacct_data_raw = pd.read_csv(StringIO(sacct_output), sep='|')

        sacct_data_raw.set_index('JobID', inplace=True)

        sacct_data = pd.DataFrame()
        for job_id in job_ids:
            if job_id in sacct_data_raw.index and f'{job_id}.batch' in sacct_data_raw.index:
                sacct_data.loc[job_id,'JobName'] = sacct_data_raw.loc[job_id,'JobName']
                sacct_data.loc[job_id,'User'] = sacct_data_raw.loc[job_id,'User']
                sacct_data.loc[job_id,'State'] = sacct_data_raw.loc[job_id,'State']
                sacct_data.loc[job_id,'ExitCode'] = sacct_data_raw.loc[job_id,'ExitCode']
                sacct_data.loc[job_id,'Elapsed'] = sacct_data_raw.loc[job_id,'Elapsed']
                sacct_data.loc[job_id,'Timelimit'] = sacct_data_raw.loc[job_id,'Timelimit']
                sacct_data.loc[job_id,'NodeList'] = sacct_data_raw.loc[f'{job_id}.batch','NodeList']
                sacct_data.loc[job_id,'Partition'] = sacct_data_raw.loc[job_id,'Partition']
                sacct_data.loc[job_id,'ReqMem'] = sacct_data_raw.loc[job_id,'ReqMem']
                sacct_data.loc[job_id,'MaxRSS'] = sacct_data_raw.loc[f'{job_id}.batch','MaxRSS']
                sacct_data.loc[job_id,'AllocCPUS'] = sacct_data_raw.loc[f'{job_id}.batch','AllocCPUS']
                sacct_data.loc[job_id,'MinCPU'] = sacct_data_raw.loc[f'{job_id}.batch','MinCPU']
                sacct_data.loc[job_id,'AveCPU'] = sacct_data_raw.loc[f'{job_id}.batch','AveCPU']
                sacct_data.loc[job_id,'MaxDiskRead'] = sacct_data_raw.loc[f'{job_id}.batch','MaxDiskRead']
                sacct_data.loc[job_id,'MaxDiskWrite'] = sacct_data_raw.loc[f'{job_id}.batch','MaxDiskWrite']
                sacct_data.loc[job_id,'Start'] = sacct_data_raw.loc[f'{job_id}.batch','Start']
                sacct_data.loc[job_id,'End'] = sacct_data_raw.loc[f'{job_id}.batch','End']

            else:
                logger.info(f'sacct data not pulled for {job_id}')
        
        self.report = self.report.merge(sacct_data, left_on='JobID', right_index=True)
                

    def update_runtime_with_log_info(self):
        from datetime import datetime
        def extract_elapsed_time(row):
            try:
                start_dt = datetime.strptime(row['Start'], '%Y-%m-%d %H:%M:%S,%f')
            except:
                logger.warning(f"Unable to extract start time\n{row}")
                start_dt = datetime.strptime('1999-01-01 01:01:01,123', '%Y-%m-%d %H:%M:%S,%f')

            try:
                logger.warning(f"Unable to extract end time\n{row}")
                end_dt = datetime.strptime(row['End'], '%Y-%m-%d %H:%M:%S,%f')
            except:
                end_dt = datetime.strptime('2099-01-01 01:01:01,123', '%Y-%m-%d %H:%M:%S,%f')
            delta = end_dt - start_dt
            hours = delta.seconds // (3600)
            minutes = (delta.seconds - hours*3600) // 60
            seconds = (delta.seconds - hours*3600 - minutes*60) 
            return f"{hours}:{minutes}:{seconds}"
        self.report['Elapsed'] = self.report.apply(extract_elapsed_time, axis = 1)

    def read_logfile(self):
        if self.blob is not None:
            print("Reading blob")
            with self.blob.open("r") as f:
                for line in f.readlines():
                    self.read_logfile_line(line)
        else:   
            for line in open(self.logfile):
                self.read_logfile_line(line)
        
        workflow_data = {key:pd.DataFrame(self.task_logs_dict[key]).T.reset_index().rename(columns={'index': 'TaskID'}) for key in self.task_logs_dict.keys()}
        for key, df in workflow_data.items():
            df['WorkflowID'] = self.workflow_id
            df['SubworkflowID'] = key
            if self.mode == 'local':
                df = df[['WorkflowID','SubworkflowID', 'TaskID', 'TaskName', 'JobID', 'Status', 'Start', 'End']]
            else:
                df = df[['WorkflowID','SubworkflowID', 'TaskID', 'TaskName', 'JobID', 'Status']]
                
            workflow_data[key] = df
        if len(workflow_data.values()) == 0:
            if self.mode == 'local':
                self.report = pd.DataFrame(columns=['WorkflowID','SubworkflowID', 'TaskID', 'TaskName', 'JobID', 'Status', 'Start', 'End'])
            else:
                self.report = pd.DataFrame(columns=['WorkflowID','SubworkflowID', 'TaskID', 'TaskName', 'JobID', 'Status'])
            return
        self.report = pd.concat(workflow_data.values())
        
        if self.workflow_name is not None:
            self.report.insert(0, 'WorkflowName', self.workflow_name)

        print(self.report.head())
        if self.mode == 'slurm':
            self.update_report_with_sacct_data()
        if self.mode == 'local':
            self.update_runtime_with_log_info()
        print(self.report.head())

    def save_report(self, odir: str):
        """Save report to specified output directory. 
        File is saved as <workflow_name>.log_report.tsv


        Args:
            odir (str): Output directory
        """
        ofile = os.path.join(odir, f'{self.workflow_name}.log_report.tsv')
        self.report.to_csv(ofile, sep='\t', index=False)

    def gantt_chart(self, ax=None):
        from matplotlib.patches import Patch
        if ax is None:
            fig, ax = plt.subplots(1, figsize=(16,6))
        self.report['Start'] = self.report['Start'].apply(lambda x: datetime.strptime(x, '%Y-%m-%dT%H:%M:%S'))
        self.report['End'] = self.report['End'].apply(lambda x: datetime.strptime(x, '%Y-%m-%dT%H:%M:%S'))
        start = self.report['Start'].min()
        self.report['rel_start_hrs'] = self.report['Start'].apply(lambda sample_start: (sample_start - start).total_seconds() / 3600)
        self.report['rel_end_hrs'] = self.report['End'].apply(lambda sample_end: (sample_end - start).total_seconds() / 3600)

        self.report['duration_hrs'] = self.report.apply(lambda row: (row.End - row.Start).total_seconds() / 3600, axis=1)


        colordict ={}
        for i, task in enumerate(self.report['TaskName'].unique()):
            colordict[task] = plt.colormaps['tab10'].colors[i]
        
        self.report['TaskColor'] = self.report['TaskName'].map(colordict)
        self.report.sort_values(by='rel_start_hrs', inplace=True, ascending=False)
        self.report.reset_index(inplace=True)

        # Make plot
        ax.barh(self.report.TaskID, self.report.duration_hrs, left=self.report.rel_start_hrs, color=self.report.TaskColor)
        # Legend
        legend_elements = [Patch(facecolor=colordict[i], label=i)  for i in colordict]
        ax.legend(handles=legend_elements)

        # Text
        for idx, row in self.report.iterrows():
            ax.text(row.rel_end_hrs+0.02, idx, 
                    f"{(row.rel_end_hrs - row.rel_start_hrs):.2f} hrs", 
                    va='center')
        ax.tick_params(left = False, right = False , labelleft = False ,)
        #Labels
        ax.set_xlabel('Time Elapsed (hours)')
        ax.set_title(f'{self.workflow_name}')
        return ax

def read_table(table_filename: str, sheet_name:str=None, header=0) -> pd.DataFrame:
    extension = os.path.splitext(table_filename)[1]
    extension = extension.lower()
    if extension in (".tsv", ".txt"):
        return pd.read_csv(table_filename, sep = '\t', header=header)
    elif extension in (".csv"):
        return pd.read_csv(table_filename, header=header)
    elif extension in (".xlsx", ".xls"):
        return pd.read_excel(table_filename, sheet_name=sheet_name, header=header)
    else:
        raise ValueError("Unknown file type. File type must be one of [.csv, .txt, .tsv, .xls, .xlsx]")


def split_gcp_filename(gcp_filename):
    split = gcp_filename.lstrip('gs://').split('/')
    bucket = split[0]
    filename = ('/').join(split[1:])
    return bucket, filename

def download_gcp_blob(bucket_name, blob_name):
    from google.cloud import storage
    """Downloads a blob into memory."""
    # The ID of your GCS bucket
    # bucket_name = "your-bucket-name"

    # The ID of your GCS object
    # blob_name = "storage-object-name"

    storage_client = storage.Client()

    bucket = storage_client.bucket(bucket_name)

    # Construct a client side representation of a blob.
    # Note `Bucket.blob` differs from `Bucket.get_blob` as it doesn't retrieve
    # any content from Google Cloud Storage. As we don't need additional data,
    # using `Bucket.blob` is preferred here.
    blob = bucket.blob(blob_name)
    return blob

def build_workflow_log_table(workflow_fof: str=None, workflow_list: list=None, workflow_root: str=None) -> pd.DataFrame:
    """Reads the workflow FOF list
    Workflow FOF optionally tab separated, contains workflow name in first column

    Args:
        workflow_log_list (str): Workflow FOF

    Returns:
        pd.DataFrame: parsed workflow files
    """
    def get_workflow_name(logfile):
        workflow_name = os.path.basename(logfile)
        workflow_name = os.path.splitext(workflow_name)[0]
        return workflow_name
    workflow_logs_from_file = pd.DataFrame(columns=['WorkflowName', 'WorkflowLog'])
    workflow_logs_from_list = pd.DataFrame(columns=['WorkflowName', 'WorkflowLog'])
    workflow_logs_from_root = pd.DataFrame(columns=['WorkflowName', 'WorkflowLog'])
    
    if workflow_fof:
        workflow_logs_from_file = read_table(workflow_fof, header=None)
        if len(workflow_logs_from_file.columns) == 2:
            workflow_logs_from_file.columns = ['WorkflowName', 'WorkflowLog']
        else:
            workflow_logs_from_file.columns = ['WorkflowLog']
            workflow_logs_from_file['WorkflowName'] = workflow_logs_from_file['WorkflowLog'].apply(get_workflow_name)
            workflow_logs_from_file = workflow_logs_from_file[['WorkflowName', 'WorkflowLog']]
    if workflow_list and workflow_list != []:
        workflow_logs_from_list = pd.DataFrame(workflow_list, columns=['WorkflowLog'])
        workflow_logs_from_list['WorkflowName'] = workflow_logs_from_list['WorkflowLog'].apply(get_workflow_name)
    if workflow_root:
        workflow_logs_from_root = build_workflow_table_from_root(workflow_root)

    workflow_logs = pd.concat([workflow_logs_from_file, workflow_logs_from_list, workflow_logs_from_root])
    workflow_logs = workflow_logs[['WorkflowName', 'WorkflowLog']]
    return workflow_logs


# def get_workflow_id(workflow_filename):
#     match = re.search(r'workflow\.(\w+)\-', workflow_filename)
#     if match:
#         return match.group(1)
#     return None

def walklevel(path, depth = 1):
    """It works just like os.walk, but you can pass it a level parameter
       that indicates how deep the recursion will go.
       If depth is 1, the current directory is listed.
       If depth is 0, nothing is returned.
       If depth is -1 (or less than 0), the full depth is walked.
       Source : https://gist.github.com/TheMatt2/faf5ca760c61a267412c46bb977718fa
    """
    # If depth is negative, just walk
    # Not using yield from for python2 compat
    # and copy dirs to keep consistant behavior for depth = -1 and depth = inf
    if depth < 0:
        for root, dirs, files in os.walk(path):
            yield root, dirs[:], files
        return
    elif depth == 0:
        return

    # path.count(os.path.sep) is safe because
    # - On Windows "\\" is never allowed in the name of a file or directory
    # - On UNIX "/" is never allowed in the name of a file or directory
    # - On MacOS a literal "/" is quitely translated to a ":" so it is still
    #   safe to count "/".
    base_depth = path.rstrip(os.path.sep).count(os.path.sep)
    for root, dirs, files in os.walk(path):
        yield root, dirs[:], files
        cur_depth = root.count(os.path.sep)
        if base_depth + depth <= cur_depth:
            del dirs[:]

def build_workflow_table_from_root(root, max_depth=2):
    log_files = []
    for path, _, files in walklevel(root, max_depth):
        for name in files:
            if  name.endswith('.log'):
                # workflow_id = get_workflow_id(name)
                relative_location = str(path).removeprefix(root).strip('/')
                # relative_location = path[len(root) + 1:]
                workflow_name = relative_location.replace('/', '.')
                # workflow_name = f'{workflow_name}.{workflow_id}' if workflow_id is not None else workflow_name
                log_files.append([workflow_name, os.path.join(path, name)])
    workflow_table = pd.DataFrame(log_files,columns=['WorkflowName', 'WorkflowLog'])
    return workflow_table


def build_workflow_reports(workflow_fof: str, workflow_list: list, workflow_root: str, mode: str, odir: str = './', name='cromreporter', save_individual=True):
    """Builds the workflow reports from the cromwell logs. Optionally queries slurm for additional
    metrics if mode=slurm. Each workflow report is saved individually. A concatenated report is also 
    generated and saved.

    Args:
        workflow_log_list (str): Filenames of workflow logs. Optionally tab separated with sample id in first column. No headers
        mode (str): Job manager. If slurm then additional metrics are generated.
        odir (str, optional): output directory. Defaults to './'.
    """
    workflow_logs = build_workflow_log_table(workflow_fof=workflow_fof, workflow_list=workflow_list, workflow_root=workflow_root)

    reports = []
    reporters = []
    current_date = datetime.now()
    date_suffix = f"{current_date.year}{current_date.month:02d}{current_date.day:02d}"
    full_report_ofile = os.path.join(odir, f'{name}.workflow_report_{date_suffix}.tsv')
    for _, row in workflow_logs.iterrows():
        logger.info(row['WorkflowName'])
        cromreporter = CromReporter(logfile=row['WorkflowLog'], workflow_name=row['WorkflowName'], mode=mode)
        reporters.append(cromreporter)

        if save_individual:
            cromreporter.save_report(odir)

        reports.append(cromreporter.report)
        cromreporter.report.to_csv(full_report_ofile, mode='a', header=not os.path.exists(full_report_ofile), index=False, sep='\t')
    full_report = pd.read_csv(full_report_ofile, sep='\t' )
    # if len(reports) > 1:
    #     logger.info('Creating final report')
    #     full_report = pd.concat(reports,ignore_index=True)

    #     # basename = os.path.basename(workflow_logs_fof)
    #     # basename = os.path.splitext(basename)[0]

    #     current_date = datetime.now()
    #     date_suffix = f"{current_date.year}{current_date.month:02d}{current_date.day:02d}"


    #     ofile = os.path.join(odir, f'{name}.workflow_report_{date_suffix}.tsv')
    #     full_report.to_csv(ofile, sep='\t', index=False)
    return full_report, reporters
    
def plot_runtime_distribution(report: pd.DataFrame,task:str, ax=None):
    report = report[report['TaskName'] == task].copy()
    if ax is None:
        fig, ax = plt.subplots(1, figsize=(16,6))
    report['Elapsed DT'] = report['Elapsed'].apply(lambda elapsed: datetime.strptime(elapsed, '%H:%M:%S') - datetime.strptime('00:00:00', '%H:%M:%S'))
    report['Elapsed (hours)'] = report['Elapsed DT'].apply(lambda td: td.total_seconds() / 3600)

    report['Elapsed (hours)'].hist(bins=20, ax=ax)
    ax.set_title(f'{task} runtime distribution')
    ax.set_xlabel('Time Elapsed (hours)')
    ax.set_ylabel('Number of runs')
    stats = []
    for i, data in report['Elapsed (hours)'].describe().round(3).items():
        stats.append(f'{i:<7} {int(data) if i =="count" else data}')
    stats = '\n'.join(stats)
    ax.text(.85,.85, stats, horizontalalignment='left',verticalalignment='center', transform=ax.transAxes,  color='gray', fontfamily='monospace')
    return ax

def plot_memory_distribution(report: pd.DataFrame,task:str, ax=None):
    report = report[report['TaskName'] == task].copy()
    def convert_to_gb(size):
        unit = size[-1]
        value = float(size[:-1])
        if unit == 'K':
            return value / (1024 * 1024)
        elif unit == 'M':
            return value / 1024
        elif unit == 'G':
            return value
        elif unit == 'T':
            return value * 1024


    report['MaxRSS (GB)'] = report['MaxRSS'].apply(convert_to_gb)
    if ax is None:
        fig, ax = plt.subplots(1, figsize=(16,6))


    report['MaxRSS (GB)'].hist(bins=20, ax=ax)
    ax.set_title(f'{task} memory usage distribution')
    ax.set_xlabel('MaxRSS (GB)')
    ax.set_ylabel('Number of runs')
    stats = []
    for i, data in report[report['TaskName'] ==task]['MaxRSS (GB)'].describe().round(3).items():
        stats.append(f'{i:<7} {int(data) if i =="count" else data}')
    stats = '\n'.join(stats)
    ax.text(.85,.85, stats, horizontalalignment='left',verticalalignment='center', transform=ax.transAxes,  color='gray', fontfamily='monospace')
    return ax

def parse_args():
    parser = argparse.ArgumentParser(description='Generates a report of all subjobs')
    parser.add_argument('-n', '--name', action='store', dest='name', default='cromreporter')
    parser.add_argument('-w', '--workflow_logs', action='store', dest='workflow_logs', nargs='*', help='workflow logs')
    parser.add_argument('-wl', '--workflow_log_list', action='store', dest='workflow_log_list', required=False,help='List of workflow logs')
    parser.add_argument('-wr', '--workflow_root',  action='store', dest='workflow_root', required=False,help='Root directory of workflow logs. Directories up to 2 levels below are searched for logs')
    parser.add_argument('-m', '--mode', choices=['slurm', 'local', ], dest='mode', default='slurm', help='Job manager. If slurm then additional metrics are generated.')
    parser.add_argument('-s', '--save-individual-reports', dest='save_individual', action=argparse.BooleanOptionalAction, default=False, help='whether to save individual workflow reports')
    parser.add_argument('-a', '--analysis', choices=['runtime', 'memory', 'gantt'], nargs='*', default=[], help='Analysis to run in addition to the tsv report')
    parser.add_argument('-odir', '--output_directory', action='store', dest='odir', default='./', help='Output directory. Directory is created if it does not exist.')
    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    os.makedirs(args.odir, exist_ok=True)

    full_report, reporters = build_workflow_reports(workflow_fof=args.workflow_log_list, workflow_list=args.workflow_logs, workflow_root = args.workflow_root, mode=args.mode, odir=args.odir, name=args.name, save_individual= args.save_individual)
    current_date = datetime.now()
    date_suffix = f"{current_date.year}{current_date.month:02d}{current_date.day:02d}"
    if 'runtime' in args.analysis and full_report is not None:
        logger.info('Creating runtime distribution plots')
        task_list = full_report['TaskName'].unique()
        runtime_filename = os.path.join(args.odir, f'{args.name}.runtime_distribution_{date_suffix}.pdf')
        with PdfPages(runtime_filename) as pdf:
            for task in task_list:
                ax = plot_runtime_distribution(full_report, task)
                pdf.savefig()  # saves the current figure into a pdf page
                plt.close()
    if 'memory' in args.analysis and full_report is not None:
        logger.info('Creating memory distribution plots')
        task_list = full_report['TaskName'].unique()
        memory_filename = os.path.join(args.odir, f'{args.name}.memory_distribution_{date_suffix}.pdf')
        with PdfPages(memory_filename) as pdf:
            for task in task_list:
                ax = plot_memory_distribution(full_report, task)
                pdf.savefig()  # saves the current figure into a pdf page
                plt.close()
    if 'gantt' in args.analysis:
        logger.info('Creating Gantt Charts')
        gantt_filename = os.path.join(args.odir, f'{args.name}.gantt_charts_{date_suffix}.pdf')
        with PdfPages(gantt_filename) as pdf:
            for crep in reporters:
                ax = crep.gantt_chart()
                pdf.savefig()  # saves the current figure into a pdf page
                plt.close()
        pass
if __name__=="__main__":
    main()

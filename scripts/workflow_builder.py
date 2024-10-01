import argparse
import subprocess
import os
import random
import string
from collections import defaultdict
import multiprocessing as mp

INPUT_TAG = "INPUT_TAG"
CROMWELL_TAG = "CROMWELL_TAG"
CONFIG_TAG = "CONFIG_TAG"
WDL_TAG = "WDL_TAG"
OPTIONS_TAG="OPTIONS_TAG"

SWARM_FILE_TAG = "SWARM_FILE_TAG"
SWARM_TIME_TAG = "TIME_TAG"
SWARM_MODULE_TAG = "MODULE_TAG"



BASE_CROMWELL_INVOCATION = "java " + CONFIG_TAG + " -jar " +  CROMWELL_TAG + " run -i " + INPUT_TAG + + " -o " + OPTIONS_TAG +" WDL_TAG"

BASE_SWARM_SUBMIT = "swarm --verbose 3 " +  SWARM_MODULE_TAG + " " + SWARM_TIME_TAG + " -t 3 -g 12 -f " + SWARM_FILE_TAG

def get_verified_absolute_path(path: str) -> str:
    """Verify and return absolute path of argument.
    Args:
        path : Relative/absolute path
    Returns:
        Absolute path
    """
    installed_path = os.path.abspath(path)
    if not os.path.exists(installed_path):
        raise RuntimeError("The requested path does not exist:{}".format(installed_path))
    return installed_path

def get_random_string(length: int) -> str:
    """
    Return a string of <length> random uppercase letters
    """
    return "".join(random.choice(string.ascii_uppercase) for i in range(0, length))

def write_to_tmp_file(filename: str, s: str) -> None:
    with open(filename, "w") as ofi:
        ofi.write(s)
    return

class WorkResult:
    def __init__(self):
        self.cromwell_log: str = None
        self.std_err_files: list = None
        self.std_out_files: list = None
        self.return_code: int = -1
        self.run_output_path: str = None
        self.inputs_json: str = "{}"
        self.outputs_json: str = "{}"

class WorkInstance:
    def __init__(self):
        self.cromwell_config_path: str = None
        self.cromwell_path: str = None
        self.wdl_path: str = None
    
        self.inputs_directory_path: str = None
        self.input_template_map: defaultdict(str) = {}
        self.input_template_json: str = None ## JSON string
        self.input_json: str = "{}" ## JSON string
        self.input_json_path: str = None
        self.rand_id: str = None
        ## A string of the form "java CONFIG_TAG -jar CROMWELL_TAG run -i INPUT_TAG WDL_TAG"
        ## should be filled by the CromRunner instance
        self.cromwell_invocation: str = None
        self.stderr: str = None
        self.stdout: str = None

    def set_rand_id(self):
        self.rand_id = get_random_string(16)

    def set_stderr_stdout(self):
        if self.rand_id is None:
            self.set_rand_id()
        p = get_verified_absolute_path(self.inputs_directory_path)
        self.stderr = p + "/" + self.rand_id + ".stderr.txt"
        self.stdout = p + "/" + self.rand_id + ".stdout.txt"
    
    def create_inputs(self):
        """
        Replace the placeholder tags in input_json with values from the input_template_map dictionary
        and write the resulting JSON to a temporary file with a unique name.
        """
        ## Get a random ID to use for writing inputs to a file.
        if self.rand_id is None:
            self.set_rand_id()
        ## Substitute inputs from a single manifest line for templated tags
        ## in inputs.
        self.input_json = self.input_template_json
        for i in self.input_template_map:
            self.input_json = self.input_json.replace( "<" + i + ">", self.input_template_map[i])
        
        p = get_verified_absolute_path(self.inputs_directory_path)
        outfile = p + "/" + self.rand_id + ".inputs.json"
        self.input_json_path = outfile
        with open(self.input_json_path, "w") as output:
            output.write(self.input_json)

        get_verified_absolute_path(self.input_json_path)
        return self.input_json_path

    def create_run_string(self):
        """
        Create a basic cromwell CLI invocation that will be passed to each WorkInstance.
        Creates the paths to cromwell, config, a tag for the input file, and a tag for the WDL
        """
        if self.cromwell_invocation is None:
            print("Error: cromwell base template must not be empty")
            exit(9)
        self.cromwell_invocation = self.cromwell_invocation.replace("INPUT_TAG", self.input_json_path)
        return self.cromwell_invocation

    def prepare_run(self):
        self.create_inputs()
        self.set_stderr_stdout()
        self.create_run_string()
        return None
    
    def run(self) -> WorkResult:
        self.prepare_run()
        try:
            with open(self.stdout, "w") as stdout_file, \
                open(self.stderr, "w") as stderr_file:
                subprocess.call(self.cromwell_invocation, shell=True, stdout=stdout_file, stderr=stderr_file)
        except(KeyboardInterrupt, Exception):
            raise KeyboardInterrupt
        return None

    def output(self):
        return

def _wrapper_func(task: WorkInstance):
    try:
        task.run()
    except(KeyboardInterrupt, Exception):
        raise KeyboardInterrupt
    return None

class CromRunner:
    """
    A server class for running multiple instances of cromwell.
    Manages a list of WorkInstances, each of which is an instantiation of Cromwell
    running a single WDL on a single input, taken from the lines of the input manifest file.
    """
    def __init__(self):
        self.available_backends = ["local", "swarm"]
        self.server_max_runtime_minutes: int = 24 * 60
        self.slurm_base_command: str = None
        
        self.swarm_base_command: str = None
        self.swarm_file: str = None
        self.swarm_submit: str = None
        self.swarm_modules: str = None

        self.work_instances: list = []
        self.max_simultaneous_instances: int = None
        self.results: list
        self.backend: str = "local"
        self.nthreads: int = 4

        self.input_manifest_path: str
        self.input_manifest_delim: str = ","
        self.input_template_path: str = None
        self.input_template: str = None
        self.wdl_path: str = None

        self.cromwell_path: str = None
        self.cromwell_config_path: str = None
        self.cromwell_invocation_template: str = BASE_CROMWELL_INVOCATION
        self.cromwell_invocation: str = ""

        self.input_tmp_dir_path: str = None
        self.dir_prefix: str = None

        self.input_manifest_header: list = None

        self.max_input_lines = 64000

    def init(self, args):
        self.cromwell_path = get_verified_absolute_path(args.cromwell_path)
        self.input_manifest_path = get_verified_absolute_path(args.manifest)
        self.nthreads = args.threads
        self.wdl_path = get_verified_absolute_path(args.wdl)
        self.input_template_path = get_verified_absolute_path(args.template)
        self.cromwell_config_path = "" if args.config is None else get_verified_absolute_path(args.config)
        self.swarm_modules_string = "" if args.modules is None else "--module " + args.modules
        self.dir_prefix = args.prefix

    def get_config_target(self):
        if self.cromwell_config_path is not None:
            return "-Dconfig.file=" + get_verified_absolute_path(self.cromwell_config_path)
        return ""

    def create_tmp_dir(self):
        rand_one = get_random_string(8)
        rand_two = get_random_string(8)
        rand_three = get_random_string(8)
        
        self.input_tmp_dir_path = "-".join([ self.dir_prefix, rand_one, rand_two, rand_three])
        os.mkdir(self.input_tmp_dir_path)
        return self.input_tmp_dir_path
        
    def create_cromwell_base_instantiation(self):
        """
        Create a basic cromwell CLI invocation that will be passed to each WorkInstance.
        Creates the paths to cromwell, config, a tag for the input file, and a tag for the WDL
        """
        self.cromwell_invocation = self.cromwell_invocation_template
        self.cromwell_invocation = self.cromwell_invocation.replace("CONFIG_TAG", self.get_config_target())
        self.cromwell_invocation = self.cromwell_invocation.replace("CROMWELL_TAG", self.cromwell_path)
        self.cromwell_invocation = self.cromwell_invocation.replace("WDL_TAG", self.wdl_path)
        return self.cromwell_invocation
    
    def load_inputs_template(self):
        with open(self.input_template_path, "r") as ifi:
            self.input_template = ifi.read()
        return self.input_template

    def create_work_instances(self):
        
        ## Create the base cromwell CLI invocation
        crom_base = self.create_cromwell_base_instantiation()

        ## Set up the input template header
        ## and iterate over each line in the input manifest,
        ## creating a new work instance for each line
        input_lines_read = 0
        with open(self.input_manifest_path, "r") as manifest:
            for line in manifest:
                line = line.strip()
                splits = line.split(self.input_manifest_delim)
                if input_lines_read == 0:
                    self.input_manifest_header = splits
                else:
                    work = WorkInstance()
                    work.inputs_directory_path = self.input_tmp_dir_path
                    work.cromwell_path = self.cromwell_path
                    work.wdl_path = self.wdl_path

                    work.input_template_json = self.load_inputs_template()

                    work.cromwell_config_path = self.get_config_target()

                    work.cromwell_invocation = crom_base

                    for i in range(0, len(self.input_manifest_header)):
                         work.input_template_map[self.input_manifest_header[i]] = splits[i]
                    self.work_instances.append(work)
                input_lines_read += 1
                if input_lines_read > self.max_input_lines:
                    print("ERROR: number of inputs exceeds maximum allowable number (" + self.max_input_lines + ").")
                    exit(9)
        print("Loaded " + str(len(self.work_instances)) + " work units.")

    def run_slurm(self):
        """
        Run a CromRunner server that submits one SLURM job
        per WorkInstance.
        """
        return

    def run_local(self):
        pool = mp.Pool(self.nthreads)
        try:
            ret = pool.map_async(_wrapper_func, self.work_instances).get(10000)
        except (KeyboardInterrupt):
            pool.terminate()
            exit(9)
        return

    def create_swarm_submit_string(self):
        ## TODO: replace time with a parameter, link to CLI, and set a reasonable default
        swarm_submit_string = BASE_SWARM_SUBMIT.replace(SWARM_TIME_TAG, "--time 36:00:00")
        swarm_submit_string = swarm_submit_string.replace(SWARM_MODULE_TAG, self.swarm_modules_string)
        return swarm_submit_string

    def run_swarm(self):
        """
        Run a CromRunner server that launches a SWARM job to handle
        all WorkInstances.
        """

        ## Step one: write the swarm file
        self.swarm_file = get_verified_absolute_path(self.input_tmp_dir_path) + "/swarm_tasks.txt"
        with open(self.swarm_file, "w") as sfi:
            for work in self.work_instances:
                work.prepare_run()
                sfi.write(work.cromwell_invocation + "\n")
        os.chmod(self.swarm_file, 0o666)

        ## Step two: write a submission script for the swarm job.
        submit_string = self.create_swarm_submit_string()
        self.swarm_submit = get_verified_absolute_path(self.input_tmp_dir_path) + "/swarm_submit.sh"
        submit_string = submit_string.replace(SWARM_FILE_TAG, self.swarm_file)
        with open(self.swarm_submit, "w") as sfi:
            sfi.write("#!/usr/bin/env bash\n")
            sfi.write(submit_string + "\n")
        os.chmod(self.swarm_submit, 0o777)
        return self.swarm_submit

    def run_gcp(self):
        """
        Launch a CromRunner server that launches jobs on Google Cloud.
        """
        return

    def __uptime(self) -> int:
        """
        Get the current uptime in minutes
        """
        return


def run_cromwell(work_instances: list) -> list:
    """
    Given a list of WorkInstances, run cromwell, once per instance.
    Return a list of WorkResults, one per instance
    """
    return

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-C", "--cromwell-path", type=str, dest="cromwell_path", default="cromwell.jar", required=False)
    parser.add_argument("-t", "--input-template", dest="template", help="Template file for WDL inputs.", required=True, type=str)
    parser.add_argument("--config", dest="config", help="An optional Cromwell config file.", required=False, type=str, default=None)
    parser.add_argument("-w", "--wdl", dest="wdl", help="A WDL script that defines a workflow.", required=True, type=str)
    parser.add_argument("-i", "--input-manifest", dest="manifest", help="A file containing inputs per instantiation, one per line, delimiter-separated, with a header that corresponds to templated input tags.", required=True, type=str)
    parser.add_argument("-d", "--delimiter", help="The delimiter used for input args in the input-file.", default=",", type=str, required=False)
    parser.add_argument("-n", "--num-concurrent-cromwells",dest="threads", help="The number of cromwell instances to spawn at one time (only applicable on local backend).", required=False, type=int, default=4)
    parser.add_argument("-B", "--backend", dest="backend", required=False, help="A backend to use for launching each cromwell instance (local or swarm). [local]")
    parser.add_argument("--modules", dest="modules", help="The modules to load when running with the SWARM backend.", type=str, default=None)
    parser.add_argument("--prefix", dest="prefix", help="A prefix to use for the CromRunner inputs directory.", type=str, default="INPUTS-CROMRUNNER")
    return parser.parse_args()

if __name__ == "__main__":

    args = get_args()
    runner = CromRunner()
    runner.init(args)
    runner.create_tmp_dir()
    invoc = runner.create_cromwell_base_instantiation()
    runner.create_work_instances()

    if args.backend == "swarm":
        runner.run_swarm()
    else:
        runner.run_local()




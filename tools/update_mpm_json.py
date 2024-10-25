import json
import argparse
import re
import glob
import os


def find_latest_checkpoint(result_dir):
    file_paths = glob.glob(f"{result_dir}/particles*.h5")
    files_names = [file_path.split("/")[-1] for file_path in file_paths]

    latest_checkpoint = None

    for file_name in files_names:
        # Use regex to find the checkpoint number in the file name
        match = re.search(r'(\d+)\.h5$', file_name)
        if match:
            checkpoint = int(match.group(1))  # Convert to int to remove leading zeros
            if latest_checkpoint is None or checkpoint > latest_checkpoint:
                latest_checkpoint = checkpoint

    return latest_checkpoint


class MPMResume:
    def __init__(self, json_path, result_dir):
        # Get checkpoint file names (particles h5 files)
        self.result_dir = result_dir
        self.file_paths = glob.glob(f"{result_dir}/particles*.h5")
        self.file_names = [file_path.split("/")[-1] for file_path in self.file_paths]

        self.file_name_info = []
        for file_name in self.file_names:
            pattern = re.compile(r'particles-(\d+)_(\d+)-(\d+)\.h5')
            match = pattern.match(file_name)
            if match:
                mpi, n_mpis, timestep = match.groups()
                self.file_name_info.append({
                    'filename': file_name,
                    'mpi': f"{mpi}_{n_mpis}",
                    'timestep': timestep
                })

        # Open and load the JSON file as a dictionary
        with open(json_path, 'r') as file:
            self.mpm_json = json.load(file)

        self.json_path = json_path

        # Get final simulation steps
        self.nsteps = self.mpm_json["analysis"]["nsteps"]

    def process_checkpoint_filename(self):
        final_checkpoint_digits = len(str(self.nsteps))

        for i, info in enumerate(self.file_name_info):
            padded_timestep = info["timestep"].zfill(final_checkpoint_digits)
            new_filename = f'particles-{info["mpi"]}-{padded_timestep}.h5'
            new_path = os.path.join(self.result_dir, new_filename)
            if self.file_paths[i] != new_path:
                os.rename(self.file_paths[i], new_path)
                print(f'Renamed: {self.file_paths[i]} -> {new_filename}')
            else:
                print(f'Skipped: {self.file_paths[i]} (has the same digits)')

    def to_latest_checkpoint(
            self, result_dir: str,
            fixed_checkpoint: int = None,
    ):
        if fixed_checkpoint is None:
            self.latest_checkpoint = find_latest_checkpoint(result_dir)
        else:
            self.latest_checkpoint = fixed_checkpoint

        # TODO: integrity check for the checkpoint detected
        # Modify the value associated with the specified key
        if self.latest_checkpoint is not None:
            self.mpm_json['analysis']['resume']['resume'] = True
            self.mpm_json['analysis']['resume']['step'] = self.latest_checkpoint
            print(f"Update mpm.json resume `true` and checkpoint to {self.latest_checkpoint}")
        else:
            print(f"Lasted checkpoint not found. mpm.json is not updated")

    def remove_particle_constraints(self):
        raise NotImplemented("This feature is not yet implemented")

    def save_updated_json(self):
        # Save the modified dictionary back to the JSON file
        with open(self.json_path, 'w') as file:
            json.dump(self.mpm_json, file, indent=4)

        print(f"Updated mpm.json saved to {self.json_path}")


if __name__ == "__main__":
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Modify a value in a JSON file.")
    parser.add_argument(
        '--json_path',
        default="/scratch1/08264/baagee/cbgeopy-scratch/simulations/sand2d-layers/mpm.json",
        type=str,
        help="The path to the JSON file.")
    parser.add_argument(
        '--result_dir',
        default="/scratch1/08264/baagee/cbgeopy-scratch/simulations/sand2d-layers/results/sand2d/",
        type=str,
        help="The directory to the mpm results to refer to.")
    parser.add_argument(
        '--update_option',
        choices=['to_latest_checkpoint', 'to_first_checkpoint', 'remove_particle_constraints'],
        default=None,
        help="The options specifying what to update")

    # Parse the arguments
    args = parser.parse_args()

    json_path = args.json_path
    result_dir = args.result_dir
    update_option = args.update_option

    # json_path = '/work2/08264/baagee/frontera/cbgeopy/examples/sand_layers-2d-random/sim-1001/mpm-resume.json'
    # result_dir = '/work2/08264/baagee/frontera/cbgeopy/examples/sand_layers-2d-random/sim-1001/results/sand2d/'
    # update_option = 'to_latest_checkpoint'

    mpm_resume = MPMResume(json_path=json_path, result_dir=result_dir)

    if update_option is not None:
        # Update particle h5 checkpoint timestep to proper digit
        mpm_resume.process_checkpoint_filename()

        if update_option == 'to_latest_checkpoint':
            mpm_resume.to_latest_checkpoint(result_dir=result_dir)
        elif update_option == 'to_first_checkpoint':
            mpm_resume.to_latest_checkpoint(result_dir=result_dir, fixed_checkpoint=0)
        elif update_option == 'remove_particle_constraints':
            mpm_resume.remove_particle_constraints()
        else:
            raise ValueError(f"Not a valid option for {update_option}")

    mpm_resume.save_updated_json()

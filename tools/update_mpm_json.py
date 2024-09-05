import json
import argparse
import re
import glob


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


class UpdateMPMjson:
    def __init__(self, json_path):
        self.latest_checkpoint = None

        # Open and load the JSON file as a dictionary
        with open(json_path, 'r') as file:
            self.mpm_json = json.load(file)

        self.json_path = json_path

    def to_latest_checkpoint(self, result_dir: str, fixed_checkpoint: int = None):
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
        default='to_latest_checkpoint',
        help="The options specifying what to update")

    # Parse the arguments
    args = parser.parse_args()

    update_mpm_json = UpdateMPMjson(json_path=args.json_path)
    if args.update_option == 'to_latest_checkpoint':
        update_mpm_json.to_latest_checkpoint(result_dir=args.result_dir)
    elif args.update_option == 'to_first_checkpoint':
        update_mpm_json.to_latest_checkpoint(result_dir=args.result_dir, fixed_checkpoint=0)
    elif args.upate_option == 'remove_particle_constraints':
        update_mpm_json.remove_particle_constraints()
    else:
        raise ValueError(f"Not a valid option for {args.update_option}")

    update_mpm_json.save_updated_json()

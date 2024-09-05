import os
import re
import glob

def remove_files_with_checkpoint_threshold(
        directory: str,
        starting_checkpoint: int,
        extensions: list,
        remove_option: str,
        checkpoint_interval=None):
    """
    This tool collects and sort all the mpm output files by the checkpoint order
         and removes the specified set of checkpoints.
         It supports removing based on the specified minimum checkpoints, and
         also supports removing based on the specified checkpoint interval.
    Args:
        directory (str): location where your mpm result files located
        starting_checkpoint (int): a minimum checkpoint to start removing
        extensions (list): a list of the output extension, e.g., ['vtu', 'pvtp', 'h5', 'vtp']
        remove_option (str): 'interval' or 'larger_or_equal'
        checkpoint_interval (int): if the remove option is 'interval', specify the interval value
    """
    for extension in extensions:
        if extension == 'pvtp':
            pattern = os.path.join(directory, f'*.{extension}')
            regex = re.compile(rf'.*?(\d+)\.{extension}$')
        elif extension == 'vtu':
            pattern = os.path.join(directory, f'*-*.{extension}')
            regex = re.compile(rf'.*-(\d+)\.{extension}$')
        else:
            pattern = os.path.join(directory, f'*-*.{extension}')
            regex = re.compile(rf'.*-(\d+)\.{extension}$')

        # Find all files matching the pattern
        files = glob.glob(pattern)

        for file in files:
            match = regex.match(file)
            if match:
                checkpoint_id = int(match.group(1))
                if remove_option == 'interval':
                    if checkpoint_interval is not None:
                        if checkpoint_id % checkpoint_interval != 0:
                            print(f"Removing {file} (interval-based)")
                            os.remove(file)
                    else:
                        raise ValueError("checkpoint_interval to remove is not specified")
                elif remove_option == 'larger_or_equal':
                    if checkpoint_id >= starting_checkpoint:
                        print(f"Removing {file} (threshold-based)")
                        os.remove(file)
                else:
                    raise ValueError(f"{remove_option} is a not a valid remove option")


if __name__ == "__main__":
    # Set your directory, checkpoint threshold, and file extension here
    directory = '/scratch1/08264/baagee/cbgeopy-scratch/simulations/fundao3d-8-3/results/sand3d/'  # Current directory
    checkpoint_threshold = 0
    extensions = ['vtu', 'pvtp', 'h5', 'vtp']
    # extensions = ['vtp']

    remove_files_with_checkpoint_threshold(
        directory, checkpoint_threshold, extensions, remove_option='interval', checkpoint_interval=10000)

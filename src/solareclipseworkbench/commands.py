import logging
import subprocess

def execute_command(command: str) -> None:
    """ Executes a command.

    Args:
        - command: Command to execute, as a string. The command should be a valid bash command.
    """

    # Execute a bash command in the shell
    logging.info(f"Executing command: {command}")

    try:
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        logging.error(f"An error occurred while executing the command: {e}")
        raise
    except Exception as e:
        logging.error(f"An unexpected error occurred: {e}")
        raise



import os
import subprocess
import sys
import logging

# Configure logging
LOG_FILE = "install_log.txt"
logging.basicConfig(
    filename=LOG_FILE,
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)


def log_and_print(message, level="info"):
    """Log a message to the log file and print it to the console."""
    print(message)
    if level == "info":
        logging.info(message)
    elif level == "error":
        logging.error(message)


def is_running_in_wsl():
    """Check if the script is running inside WSL."""
    try:
        with open("/proc/version", "r") as version_file:
            return "microsoft" in version_file.read().lower()
    except FileNotFoundError:
        return False


def ensure_wsl_installed():
    """Check if WSL is installed (only relevant outside WSL)."""
    if is_running_in_wsl():
        log_and_print("✅ Script is running inside WSL.")
        return

    try:
        log_and_print("🔄 Checking if WSL is installed...")
        subprocess.run("wsl --list", shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        log_and_print("✅ WSL is installed.")
    except subprocess.CalledProcessError:
        log_and_print("❌ WSL is not installed. Please install WSL before running this installer.", level="error")
        print("\nTo install WSL manually, run the following command in an administrator PowerShell:")
        print("\nwsl --install")
        print("\nFor more information, visit: https://aka.ms/wslinstall\n")
        input("Press any key to close this window after reading the instructions...")
        sys.exit(1)


def ensure_ubuntu_installed():
    """Ensure Ubuntu is installed (only relevant outside WSL)."""
    if is_running_in_wsl():
        log_and_print("✅ Script is running inside WSL. Skipping Ubuntu checks.")
        return

    try:
        log_and_print("🔄 Checking if Ubuntu is installed in WSL...")
        result = subprocess.run("wsl --list", shell=True, text=True, capture_output=True)
        if "Ubuntu" not in result.stdout:
            log_and_print("⚠️ Ubuntu is not installed in WSL. Installing Ubuntu...")
            subprocess.run("wsl --install -d Ubuntu", shell=True, check=True)
        log_and_print("✅ Ubuntu is installed.")
    except subprocess.CalledProcessError as e:
        log_and_print(f"❌ Error ensuring WSL setup: {e.stderr}", level="error")
        sys.exit(1)

    try:
        log_and_print("🔄 Setting Ubuntu as the default WSL distribution...")
        subprocess.run("wsl --set-default Ubuntu", shell=True, check=True)
        log_and_print("✅ Ubuntu is now the default WSL distribution.")
    except subprocess.CalledProcessError as e:
        log_and_print(f"❌ Error setting default WSL distribution: {e.stderr}", level="error")
        sys.exit(1)


def install_dependencies():
    """Install dependencies interactively in WSL."""
    ensure_wsl_installed()
    ensure_ubuntu_installed()

    # Check for required tools
    tools_to_check = {
        "hisat2": "hisat2 --version",
        "samtools": "samtools --version",
        "featureCounts": "featureCounts --version",
        "pip3": "pip3 --version"
    }

    # Determine which tools need to be installed
    tools_to_install = []
    for tool, check_cmd in tools_to_check.items():
        if not check_tool_installed(tool, check_cmd):
            tools_to_install.append(tool)

    # If all tools are installed, exit
    if not tools_to_install:
        log_and_print("✅ All dependencies are already installed!")
        return

    log_and_print("⚠️ Some dependencies are missing. Installing them in WSL...")

    # Run installation commands interactively
    run_interactive_command("sudo apt-get update -y", "Updating package lists")
    run_interactive_command(
        "sudo apt-get install -y hisat2 samtools subread python3 python3-pip",
        "Installing missing tools",
    )
    run_interactive_command(
        "pip3 install PySide6 paramiko",
        "Installing Python libraries",
    )


def check_tool_installed(tool_name, check_command):
    """Check if a tool is installed in WSL."""
    try:
        subprocess.run(check_command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        log_and_print(f"✅ {tool_name} is already installed.")
        return True
    except subprocess.CalledProcessError:
        log_and_print(f"⚠️ {tool_name} is not installed.")
        return False


def run_interactive_command(command, description):
    """Run a shell command interactively."""
    try:
        log_and_print(f"🔄 {description}")
        result = subprocess.run(command, shell=True, check=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        log_and_print(result.stdout)
        return result
    except subprocess.CalledProcessError as e:
        log_and_print(f"❌ Error during {description}: {e.stderr}", level="error")
        sys.exit(1)


if __name__ == "__main__":
    try:
        install_dependencies()
        log_and_print("✅ All dependencies are installed!")
        log_and_print("🎉 Installation completed successfully.")
    except Exception as e:
        log_and_print(f"❌ An unexpected error occurred: {e}", level="error")
        input("Press any key to close this window after reviewing the error...")

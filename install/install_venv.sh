
# get env name from var 1 (required), and path from var 2 (optional)

# Usage:
# source install_venv.sh ENV_NAME [PATH]

# Example:
# source install_venv.sh TEST_ENVIRONMENT /isiseqruns/jfreeman_tmp_home/bin/miniconda3/envs

# If no path is provided, the default path is used:
# /isiseqruns/jfreeman_tmp_home/bin/miniconda3/envs

# If the environment already exists, it will be overwritten

# check that conda is installed:
if ! command -v conda &> /dev/null
then
    echo "conda could not be found"
    echo "Please install conda and try again"
    exit 1
fi

#check that there is at least one argument:
if [ $# -eq 0 ]
  then
    echo "Virtual Environment name is required"
    echo "Usage: source install_venv.sh ENV_NAME [PATH]"
    exit 1
fi

if [ "$(basename "$PWD")" != "install" ]
then
  echo "Please cd into the install directory"
fi

# Set the name of the environment equal to the first argument
name=$1

# If a second argument is provided, set the path equal to it
if [ $# -eq 2 ]
  then
    conda_envs_path=$2
else
    # deduce the path from the conda executable path:
    conda_envs_path=$(dirname $(dirname $(which conda)))/envs
    # make the directory if it doesn't exist:
    mkdir -p $conda_envs_path
fi

# Construct the full path
full_path="${conda_envs_path}/${name}"

# Replace PREFIX_PLACEHOLDER with the desired path in the .yml file
sed "s|PREFIX_PLACEHOLDER|$full_path|" environment_template.yml > environment.yml

# Replace NAME_PLACEHOLDER with the desired name in the generated environment.yml file
sed -i "s|NAME_PLACEHOLDER|$name|" environment.yml

conda env create -f environment.yml

conda activate $name

R -e "devtools::install_github('immunogenomics/harmony')"


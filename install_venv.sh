. $HOME/miniconda/etc/profile.d/conda.sh

CONDA_PATH=$HOME/miniconda/bin/conda

if ! command -v $CONDA_PATH &> /dev/null
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
    return 1
fi

# Set the name of the environment equal to the first argument
name=$1

# If a second argument is provided, set the path equal to it
if [ $# -eq 2 ]
  then
    conda_envs_path=$2
else
    # deduce the path from the conda executable path:
    conda_envs_path=$(dirname $(dirname $CONDA_PATH))/envs
    # make the directory if it doesn't exist:
    mkdir -p $conda_envs_path
fi

# Construct the full path
full_path="${conda_envs_path}/${name}"

# Replace PREFIX_PLACEHOLDER with the desired path in the .yml file
sed "s|PREFIX_PLACEHOLDER|$full_path|" environment_template.yml > environment.yml

# Replace NAME_PLACEHOLDER with the desired name in the generated environment.yml file
sed -i "s|NAME_PLACEHOLDER|$name|" environment.yml

# First, remove "/envs" from the end of conda_envs_path
conda_bin_path="${conda_envs_path%/envs}"

# Now add "/condabin/conda" to the end
conda_bin_path="${conda_bin_path}/condabin/conda"

# Use this new variable in your sed command:

# Full path to the R executable
r_path="${full_path}/bin/R"


$CONDA_PATH env create -f environment.yml




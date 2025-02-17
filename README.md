# splicai utils
This repository contains files required to run SpliceAI genome-wide in HPC environments. The following files are included:

1. ``runSpliceAi.perl``: Overall pipeline script, invokes the others
2. ``pipeable_spliceai_to_wiggle_twoBit``: Runs spliceai for a given section of an assembly
3. ``para-nf``: A minimalist bash file to digest joblist files using nextflow.
4. ``nf-joblist.nf``: Joblist required by para-nf

# Prequisites

Before you can run these script, the following changes / installs are required:

Create python virtual environment (venv) and install the dependencies outlined in requirements.txt. Update the shebang of the ``pipeable_spliceai_to_wiggle_twoBit`` script to point to the interpreter of the python venv. 

In order to run SpliceAI, you'll need the models from the original publication (-> https://github.com/Illumina/SpliceAI/tree/master/spliceai/models). Download the models and adjust the path on line 98 in ``pipeable_spliceai_to_wiggle_twoBit`` to point towards the location of the model .h5 files on your system.

On its last line, `para-nf` references the file location of `nf-joblist.nf` directly. Update the file path for your system.

Finally, ``pipeable_spliceai_to_wiggle_twoBit`` and ``para-nf`` need to be in your PATH. Either put the into the corresponding /bin directories or symlink them there. Make sure that the execute permission is set.

# splicai utils
This repository contains files required to run SpliceAI genome-wide in HPC environments. The following files are included:

1. ``runSpliceAi.perl``: Overall pipeline script, invokes the others
2. ``pipeable_spliceai_to_wiggle_twoBit``: Runs spliceai for a given section of an assembly
3. ``para-nf``: A minimalist bash file to digest joblist files using nextflow. Used internally for upgrading to nextflow.
4. ``nf-joblist.nf``: Joblist required by para-nf

# Prequisites

Before you can run these script, the following changes / installs are required:

Install python 3.9.2 -> newer versions may work as well but have not been tested. Create python virtual environment (venv or similar) and install the dependencies outlined in requirements.txt. These package versions work with the most recent 1.3.1 version of SpliceAI. Update the shebang of the ``pipeable_spliceai_to_wiggle_twoBit`` script to point to the interpreter of the python venv.

In order to run SpliceAI, you'll need the models from the original publication (-> https://github.com/Illumina/SpliceAI/tree/master/spliceai/models). Download the models and adjust the path on line 98 in ``pipeable_spliceai_to_wiggle_twoBit`` to point towards the location of the model .h5 files on your system. Note that the licence of the models differs from the licence for the spliceai package.

On its last line, `para-nf` references the file location of `nf-joblist.nf` directly. Update the file path for your system. While SpliceAI jobs generally run very fast, the generation of genome-wide bigWig files requires time and memory (~200GB). If your HPC environment places restrictions on runtimes or memory usage, add the `--queue` parameter to the `para-nf` calls in the perl script to tell nextflow to push the job to the correct queue / partition.

Finally, ``pipeable_spliceai_to_wiggle_twoBit`` and ``para-nf`` need to be in your PATH. Either put the into the corresponding /bin directories or symlink them there. Make sure that the execute permission is set.

---
The contents of this repository were created at the hillerlab (https://github.com/hillerlab)
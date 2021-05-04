The genome scale model for rhodococcus opacus was created using carveme. For more details about carveme refer to the following paper: https://pubmed.ncbi.nlm.nih.gov/30192979/

Model generation was performed on a terminal in Jprime.Carveme is already installed in the biodesign kernel on Jprime.

To activate the kernel on terminal, type the following command: "pyenv activate biod_3.7"
If you need to reinstall carveme somewhere else look into the following documentation: http://carveme.readthedocs.io/
Reader should keep in mind that carveme needs cplex and diamond in order for it to create model.

Installing Diamond is specific for each user on Jprime, so make sure you install diamond locally on your user account.
For Diamond installing instructions look into the following: http://www.diamondsearch.org

After installing all dependencies, initiate carveme with the following command: "carveme_init"
To create R.opacus model, the following genome data was used: https://www.ncbi.nlm.nih.gov/assembly/GCF_000234335.1
The following command creates R.opacus GEM: "carve --refseq GCF_000234335.1 -u grampos -o Ropacus_carveme.xml"

After the model is created, perform gapfilling for the model using the following command: "gapfill Ropacus_carveme.xml -m M9,LB -o new_model.xml"
The model should be ready and was curated using the notebook given in this repo.

/* SJS. 
HyPhy implementation of Rate4Site. 
Currently uses JukesCantor model, in JC_aa.mdl. Matrix and frequency vector names are hard-coded, because I can't figure out how not to do that in HyPhy at the moment. 
*/


LIKELIHOOD_FUNCTION_OUTPUT = 1;
RANDOM_STARTING_PERTURBATIONS = 1;
OPTIMIZATION_PRECSION = 0.0000000001;
ACCEPT_BRANCH_LENGTHS=1;

global r;
#include "setup_helper.ibf"; // code to read setup

// Setup is read from default from file "setup.txt".
// If this file doesn't exist, specify setup file from
// command line like this:
//     HYPHYMP FEL.bf <<< "mysetup.txt"
setup = readSetup();
infile = setup["INFILE"];
outfile = setup["OUTFILE"];

#include "JC_aa.mdl"; // Can't find a way to have this file name be a variable in setup. Unfortunate. Either way, it currently also contains a frequency vector (equal) called JC_freqs.


DataSet raw_data  = ReadDataFile(infile);
DataSetFilter filtered_data = CreateFilter(raw_data, 1);

fprintf(stdout, "Step 1: Global optimization of branch lengths.\n");
Model JCfull = (JC69_t, JC_freqs, 1);
UseModel(USE_NO_MODEL);
UseModel(JCfull);
Tree full_tree = DATAFILE_TREE;
LikelihoodFunction full_ln_likfn = (filtered_data, full_tree);
Optimize(full_res, full_ln_likfn);




fprintf(stdout, "Step 2: Site-specific rate scalar computation.\n");
header = "site\trate\tlnL\n";
fprintf(outfile, CLEAR_FILE, header);

nsites = filtered_data.sites;
for (global site_count = 0; site_count < nsites; site_count = site_count+1)
{
       
    // single amino acid
    filter_string = "";
    filter_string = filter_string + (site_count); // + "-" + (site_count*3+2);
    DataSetFilter site_filter = CreateFilter(filtered_data, 1, filter_string, "", "");

    Model JCsite = (JC69_rt, JC_freqs, 1);
    UseModel(USE_NO_MODEL);
    UseModel(JCsite);
    Tree site_tree = DATAFILE_TREE;
    LikelihoodFunction site_ln_likfn = (site_filter, site_tree);

    // Two underscores on right side of constraint mean "take the value of the variable,
    // not the variable itself"
    ReplicateConstraint("this1.?.t:=this2.?.t__", site_tree, full_tree);
    Optimize(site_res, site_ln_likfn);

    site_rate = site_res[0][0];
    site_lnlk = site_res[1][0];
    line = Format(site_count+1, 0, 0) + "\t" + Format(site_rate, 5, 10) + "\t" + Format(site_lnlk, 5, 10) + "\n";
    fprintf(stdout, line);
    fprintf(outfile, line);
}

import sys
import os 

tree_file = sys.argv[1]
orthogroup = tree_file.split("/")[4].split("_")[0]

for filename in os.listdir(sys.argv[2]):

	if orthogroup in filename:
		seq_file = os.path.join(sys.argv[2], filename)
		break

with open(("paml_" + orthogroup + "_codeml.ctl"), 'w') as control_file:

	control_file.write("seqfile = " + seq_file + "\n")
	control_file.write("outfile = " + orthogroup + "_paml.out\n")
	control_file.write("treefile = " + tree_file + "\n")
	control_file.write("\n")
	control_file.write("noisy = 9\n")
	control_file.write("verbose = 0\n")
	control_file.write("runmode = -2\n")
	control_file.write("\n")
	control_file.write("\n")
	control_file.write("seqtype = 1\n")
	control_file.write("CodonFreq = 2\n")
	control_file.write("* ndata = 10\n")
	control_file.write("clock = 0\n")
	control_file.write("\n")
	control_file.write("aaDist = 0\n")
	control_file.write("\n")
	control_file.write("model = 0\n")
	control_file.write("\n")
	control_file.write("NSsites = 0\n")
	control_file.write("\n")
	control_file.write("icode = 0\n")
	control_file.write("\n")
	control_file.write("Mgene = 0\n")
	control_file.write("\n")
	control_file.write("fix_kappa = 0\n")
	control_file.write("kappa = 2\n")
	control_file.write("fix_omega = 0\n")
	control_file.write("omega = 0.4\n")
	control_file.write("\n")
	control_file.write("fix_alpha = 1\n")
	control_file.write("alpha = 0\n")
	control_file.write("Malpha = 0\n")
	control_file.write("ncatG = 3\n")
	control_file.write("\n")
	control_file.write("fix_rho = 1\n")
	control_file.write("rho = 0\n")
	control_file.write("\n")
	control_file.write("getSE = 0\n")
	control_file.write("RateAncestor = 0\n")
	control_file.write("\n")
	control_file.write("Small_Dif = .5e-6\n")
	control_file.write("*cleandata = 0\n")
	control_file.write("*fix_blength = 0\n")
	control_file.write("method = 0\n")
	

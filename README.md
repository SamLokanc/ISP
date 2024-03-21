# ISP
This python script will process the raw xls output of netMHCpan and netMHCIIpan. 

# Use:
python3 ISP.py <path/to/netMHCpan_out_file> <path/to/netMHCIIpan_out_file> <outfile.csv> <num_top_binders_ranked> <amino_acids_containing_core_mutation>

# Example:
Say you predicted the binding of a peptide containing a specific SNV : EEEEEAEEEEE (where A indicates the SNV). We'll call the netMHCpan I and II outputs netMHCpan_out.xls and netMHCIIpan_out.xls. To retrieve the top 10 binding epitopes for each MHC supertype use the following command:
python3 ISP.py netMHCpan_out.xls netMHCIIpan_out.xls out.csv 10 EAE

it is important in the <amino_acids_containing_core_mutation> to include the immediate flanking peptides as well.


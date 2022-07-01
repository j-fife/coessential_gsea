from config import module_paths
import glob
import pandas as pd
import gseapy as gp

class CoessentialGSEA:

    def __init__(self, input_file):
        self.input_file = input_file
        self.module_path_dir = module_paths[0.2]
        self.out_path = "./CoessentialGSEA_output/"
        self.modules = pd.read_csv(module_paths[0.2], index_col = "Cluster")

    def combine_results(self):
    	main_df = None
    	subdirs = glob.glob(self.out_path + "*")
    	for f in subdirs:
    		csvs = glob.glob(f + "/*.csv")
    		if len(csvs) == 1:
    			c = pd.read_csv(csvs[0])
    			if main_df is None:
    				main_df = c
    			else:
    				main_df = pd.concat([main_df, c])
    	return main_df

    def run(
        self,
        d = 0.2,
        permutation_num = 10000,
        output_name = None,
        processes = 1,
        min_size = 4,
        seed=0
    ):
        if output_name == None:
            output_name = "output_coessentialGSEA.csv"

        if input_file[-4:] == ".csv":
            rnk = pd.read_csv(self.input_file, skiprows = 1, names = [0,1])
        else:
            rnk = pd.read_csv(self.input_file, delim_whitespace = True, names = [0,1])
        rnk = rnk.dropna()
        rnk = rnk.sort_values(by = 1, ascending = False)
        if d in module_paths.keys():
            self.module_path_dir = module_paths[d]
            self.modules = pd.read_csv(module_paths[d], index_col = "Cluster")

        gene_members = set(rnk[0].values) ## creates set out of values in column index 1 of rnk
        counter = 0
        gene_sets = {}
        ## for every module, add list of their member genes to a set
        for index, row in self.modules.iterrows():
            members = row["Members"].split(" ")
            members = list(filter(lambda x : x in gene_members, members))
            gene_sets[f"module_{index}"] = members

        pre_res = gp.prerank(
            rnk=rnk,
            gene_sets = gene_sets,
            processes=processes, ## 1 process; upping number would be parallel processes
            permutation_num=permutation_num, ## reduce number to speed up testing (1000 or less) ; take score and do 10,000 permutations (like sample size)
            outdir= self.out_path,
            format='png',
            min_size = min_size,
            seed=seed
        )
        out_df = pd.read_csv(
            self.out_path + "gseapy.prerank.gene_sets.report.csv",
            index_col = "Term"
        )
        self.modules["Term"] = list(
            map(lambda x : "module_" + str(x), list(self.modules.index))
        )
        self.modules.set_index("Term", inplace = True)
        GO_columns = [
            "Top GO Terms",
            "Top GO Term p-values",
            "Top GO FDRs",
            "Top GO Term Fold Enrichments"
        ]
        out_df = out_df.join(self.modules[GO_columns])
        out_df.to_csv(output_name)


if __name__ == "__main__":
    input_file = "../data/example_data.csv"
    c = CoessentialGSEA(input_file)
    c.run(d = 0.2, permutation_num = 10)

import pandas as pd
import sys

# Load the data
input_tsv = sys.argv[1]
output_file = sys.argv[2]

def get_genotyp_numbers(blastfmt6_path:str):
    import csv
    #parsing tsv file
    with open(blastfmt6_path) as file:
        tsv_file = csv.reader(file, delimiter="\t")
        # in blast fmt 6 the 2nd entry is the database hit
        # The sequences in the database are named in this way: [GenotypNumber]_[SequenceName]
        # so spliting by "_" and taking only the first substring results in the Genotyp number
        #The number for the six relevent Segments (all but HA and NA) are added to a list
        genotype_numbers = []
        for line in tsv_file:
            genotype_numbers.append(line[1].split("_")[0])
    return genotype_numbers

def get_genotyp_string(genotyp_numbers:list):
    #turn genotyp_list into string seperated by |
    genotype ="|"
    genotype = genotype.join(genotyp_numbers)
    return(genotype)

def get_genotyp(genotyp_string:str):
    #initalize dictionary with genotyp strings as keys and Genotyps as Value
    genotypes = {
        '20|20|20|20|20|20|20|20':'A',
        '35|28|1|20|1|30|20|20':'B',
        '1|1|1|20|1|1|20|1':'C',
        '20|20|20|20|20|27|20|20':'D',
        '20|20|29|20|20|27|20|20':'E',
        '27|20|29|20|37|27|20|20':'F',
        '31|31|32|20|12|17|20|1':'G',
        '20|31|29|20|27|63|20|20':'H',
        '25|14|25|20|26|13|20|28':'I',
        '31|23|5|20|26|2|20|28':'J',
        '31|14|24|20|26|20|20|29':'K',
        '25|1|1|20|1|13|20|1':'L',
        '20|20|14|20|20|20|20|20':'M',
        '12|33|20|20|36|64|20|20':'N',
        '20|29|14|20|20|20|20|20':'O',
        '31|23|3|20|26|2|20|28':'P',
        '20|20|20|20|26|20|20|20':'Q',
        '20|20|29|20|15|27|20|20':'R',
        '20|20|29|20|20|20|20|20':'S',
        '25|1|1|20|1|27|20|20':'T',
        '20|20|29|20|15|20|20|20':'U',
        '34|1|1|20|1|1|20|1':'V',
        '20|9|14|20|16|61|20|20':'X',
        '1|1|1|20|26|1|20|1':'AA',
        '31|1|3|20|38|1|20|1':'AB',
        '4|1|1|20|1|1|20|1':'AC',
        '4|27|1|20|31|1|20|1':'AD',
        '4|1|1|20|15|1|20|1':'AE',
        '12|6|1|20|50|1|20|29':'AF',
        '12|1|1|20|50|1|20|1':'AG',
        '12|1|1|20|1|1|20|1':'AH',
        '7|1|8|20|37|1|20|29':'AI',
        '10|1|1|20|1|1|20|1':'AJ',
        '10|1|3|20|38|1|20|1':'AK',
        '1|1|1|20|11|1|20|1':'AL',
        '1|1|1|20|12|1|20|1':'AM',
        '31|1|1|20|1|1|20|1':'AN',
        '13|1|1|20|38|1|20|1':'AO',
        '45|31|32|20|32|1|20|28':'AP',
        '31|1|8|20|37|1|20|27':'AQ',
        '1|18|14|20|16|62|16|6':'AR',
        '1|1|3|20|38|1|20|1':'AS',
        '12|1|1|20|37|1|20|1':'AT',
        '19|1|14|20|21|22|20|1':'AU',
        '1|1|3|20|1|1|20|1':'AV',
        '4|23|1|20|11|1|20|27':'AW',
        '31|31|1|20|38|1|20|1':'AX',
        '1|1|24|20|38|1|20|1':'AY',
        '1|31|1|20|12|1|20|1':'AZ',
        '13|1|1|20|1|1|20|1':'BA',
        '31|1|43|20|43|1|20|43':'BB',
        '45|1|43|20|37|1|20|27':'BC',
        '31|31|8|20|26|1|20|1':'BD',
        '4|1|14|20|26|1|20|1':'BE',
        '13|1|1|20|11|1|20|1':'BF',
        '31|1|3|20|16|1|16|1':'CA',
        '32|1|3|20|38|1|20|1':'CB',
        '12|1|8|20|38|1|20|1':'CC',
        '31|1|3|20|37|1|20|29':'CD',
        '44|1|3|20|38|1|20|1':'CE',
        '31|1|3|20|15|1|20|1':'CF',
        '4|1|3|20|26|1|20|1':'CG',
        '31|1|3|20|26|1|20|1':'CH',
        '10|1|12|20|38|1|20|1':'CI',
        '46|6|32|20|26|22|47|1':'CJ',
        '4|1|3|20|38|1|20|1':'CK',
        '31|1|48|20|38|1|20|1':'CL',
        '31|49|3|20|38|1|20|1':'CM',
        '31|31|3|20|26|1|20|1':'CN',
        '31|31|3|20|38|1|20|1':'CP',
        '31|1|3|20|37|1|20|1':'CQ',
        '31|1|3|20|38|1|20|28':'CR',
        '34|1|3|20|38|1|20|51':'CS',
        '31|1|43|20|43|1|20|1':'CT',
        'X|29|22|20|37|1|20|27':'CU',
        '4|52|3|20|38|1|20|27':'DA',
        '31|53|3|20|38|1|20|1':'DB',
        '41|53|32|20|36|1|20|1':'DC',
        '41|1|3|20|26|1|20|56':'DD',
        '41|53|3|20|55|1|20|56':'DE',
        '57|1|14|20|38|1|20|1':'DG',
        '4|52|3|20|11|1|20|27':'DH',
        '14|14|3|20|38|1|20|27':'DI',
        '31|52|3|20|59|1|20|58':'DJ',
        '4|1|12|20|26|1|20|1':'DK',
        '31|1|1|20|38|1|20|1':'DL',
        '31|52|3|20|11|1|20|27':'DM',
        '4|60|3|20|26|1|20|27':'DN',

    }
    if  genotyp_string in genotypes:
        return(genotypes[genotyp_string])
    else:
        return("Unknown")

def get_genotype_from_blastfmt6(blastfmt6_path:str):
        genotyp_list = get_genotyp_numbers(blastfmt6_path)
        print(genotyp_list)
        genotyp_string = get_genotyp_string(genotyp_list)
        print(genotyp_string)
        genotype = get_genotyp(genotyp_string)
        
        return genotype


genotyp = get_genotype_from_blastfmt6(input_tsv)

print(genotyp,end="")   

with open(output_file, 'w') as f:
    f.write(genotyp)
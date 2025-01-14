#!/usr/bin/env python3

import sys
import os
import json
import time
import random
import zipfile
import argparse
import threading
import subprocess
from tqdm import tqdm
from Bio import Entrez
from pathlib import Path
import concurrent.futures
from multiprocessing import Pool

def download_genome(species, output_dir, max_retries=3, file_types='genome,gtf'):
    """
    Download genome and GTF files for a given species using NCBI datasets tool
    
    Args:
        species (str): Species name or taxon ID
        output_dir (Path): Directory to save downloaded files
        file_types: list of file types to download e.g. genome,gtf,etc
        
    """
    # Add random delay before starting download
    delay = random.uniform(1, 5)
    time.sleep(delay)
    print(f"\nProcessing {species} (after {delay:.1f}s delay)...")
    
    # Create species-specific directory
    species_dir = output_dir / species.replace(" ", "_")
    species_dir.mkdir(parents=True, exist_ok=True)
    
    # Check if the file(s) already exist
    dataset_catalog = species_dir / "ncbi_dataset" "/data/dataset_catalog.json"
    if (dataset_catalog).exists():
        with open(dataset_catalog, 'r') as json_file:
            species_accession = json.load(json_file)
            files_index = None
            for index, element in enumerate(species_accession['assemblies']):
                try:
                    list(element.keys()).index('accession')
                    files_index = index
                    break
                except ValueError:
                    continue
            
            if files_index is None:
                print('No Accession is found in:', files_index, species_accession['assemblies'], 'Trying to download again.')  
                #return (species, False)
            else:
                fasta_file = ''
                gtf_file = ''
                fasta = False
                gtf = False
                
                for file in species_accession['assemblies'][files_index]['files']:
                    file_type = file['fileType']
                    file_path = file['filePath']
                    if file_type == 'GENOMIC_NUCLEOTIDE_FASTA':
                        fasta_file = species_dir / "ncbi_dataset" "/data" / file_path
                    elif file_type == 'GTF':
                        gtf_file = species_dir / "ncbi_dataset" "/data" / file_path
                        
                if fasta_file:
                    if fasta_file.exists():
                        fasta = True
                
                if gtf_file:
                    if gtf_file.exists():
                        gtf = True
                else: #indicating no annotation is available for this species (GTF is missin in the dataset report)
                    print(f"No GTF files is avilable for: {species}")
                    gtf = True
                
                if fasta and gtf:
                    print(f'Files already exist for species: {species}\n{fasta_file} {gtf_file}')
                    return (species, True)
                else:
                    print("Files don't exist. Trying to download")
    
    #Dowload the files
    try:
        for attempt in range(max_retries):
            print(f'Attempt: {attempt}. Downloading: {species} to {species_dir}')
            # Download genome data using datasets command
            cmd = f"datasets download genome taxon \"{species}\" --include {file_types} --reference"
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True, cwd=species_dir)
            
            if result.returncode == 0:
                #Unzip downloaded folder that contains the files
                zip_path =  species_dir / "ncbi_dataset.zip"
                extract_path = species_dir / "ncbi_dataset"
                try:
                    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
                        zip_ref.extractall(extract_path)
                        # Check if download was successful by looking for the data report
                        if dataset_catalog.exists():
                            # Clean up zip file
                            zip_file.unlink()
                            print(f"Successfully downloaded data for {species}")
                            return (species, True)
                except zipfile.BadZipFile as e:
                    print(f'Error in unzipping {speceis}. Error: {e}')
                except FileNotFoundError as e:
                    print(f'Error in unzipping {speceis}. Error: {e}')
    except Exception as e:
        print(f"Error processing {species}")
    species_dir.rmdir()
    return (species, False)

def get_species_names(search_term, email="husen.umer@kaust.edu.sa", api_key=None, max_results=100, accepted_divisions=['all']):
    """
    Fetches a list of species names from NCBI taxonomy database based on a search term.
    
    Args:
        search_term (str): The term to search for in the NCBI taxonomy database.
        email (str): Your email address (required by NCBI's Entrez).
        api_key (str): NCBI API key (optional).
        max_results (int): Maximum number of results to fetch.
        
    Returns:
        list: A list of species names.
    """
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

    try:
        # Search the taxonomy database
        handle = Entrez.esearch(db="taxonomy", term=search_term, retmax=max_results)
        record = Entrez.read(handle)
        handle.close()
        
        # Get the taxonomy IDs
        tax_ids = record["IdList"]
        if not tax_ids:
            print("No species found for the given search term.")
            return []
        
        # Fetch details for each taxonomy ID
        handle = Entrez.efetch(db="taxonomy", id=",".join(tax_ids), retmode="xml")
        records = Entrez.read(handle)
        handle.close()
        
        # Extract species names
        species_names = []
        for record in records:
            try:
                if record['Division'] in accepted_divisions or accepted_divisions==['all']:
                    species_names.append(record["ScientificName"])
            except KeyError:
                pass
        return species_names

    except Exception as e:
        print(f"An error occurred: {e}")
        return []

def fetch_ncbi_species(email, max_species, accepted_divisions, species_list_file):
    """
    Fetch species names from NCBI Taxonomy database.
    
    Args:
        email (str): Your email address (required by NCBI)
        batch_size (int): Number of records to fetch per batch
        sleep_interval (float): Time to wait between requests in seconds
        species_list_file: file name to store list of species from NCBI
    
    Returns:
        list: List of species names
    """
    batch_size=10000
    sleep_interval=1
    Entrez.email = email  # Required by NCBI
    species_list = []
    divisions = set()
    # First, get the total count of species
    handle = Entrez.esearch(db="taxonomy", term="species[RANK]")
    record = Entrez.read(handle)
    total_count = int(record["Count"])
    handle.close()
    
    # Fetch species in batches
    count_results = 0
    with open(species_list_file, 'w') as species_list_outfile:
        for start in range(0, total_count, batch_size):
            try:
                # Fetch a batch of species
                handle = Entrez.esearch(db="taxonomy", 
                                    term="species[RANK]", 
                                    retstart=start, 
                                    retmax=batch_size)
                record = Entrez.read(handle)
                handle.close()
                
                # Get details for each species ID
                id_list = record["IdList"]
                handle = Entrez.efetch(db="taxonomy", id=id_list, retmode="xml")
                records = Entrez.read(handle)
                handle.close()
                
                # Extract species names
                for record in records:
                    species_name = record.get("ScientificName", "")
                    if record['Division'] in accepted_divisions or accepted_divisions==['all']:
                        if species_name:
                            species_list.append(species_name)
                            species_list_outfile.write(species_name + '\n')
                    if len(species_list) >= max_species:
                        print(f'Already obtained {max_species} species. Exiting')
                        return species_list
                print(f"Processed {len(species_list)} species so far...")
                time.sleep(sleep_interval)  # Be nice to NCBI servers
            except Exception as e:
                print(f"Error occurred at offset {start}: {str(e)}")
                continue
    
    return species_list

def main():
    parser = argparse.ArgumentParser(description="Download genome and GTF files from NCBI for multiple species")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-s", "--species", nargs="+", 
                        help="List of species names or taxon IDs")
    group.add_argument("-f", "--file", type=str, 
                        help="File containing species names or taxon IDs (one per line)")
    group.add_argument("-st", "--search_term", type=str,
                        help="Search term to identify species name and if it exists")
    group.add_argument("-l", "--list_all_species", action='store_true', default=False,
                        help="List all species available in NCBI")
    group.add_argument("-ad", "--accepted_divisions", nargs="+", default=['all'],
                        help="""List of accepted divisions including: 'Bacteria', 'Vertebrates', 'Invertebrates', 'Viruses', 'Environmental samples', 
                            'Mammals', 'Synthetic and Chimeric', 'Unassigned', 'Phages', 'Plants and Fungi', 'Primates', 'Rodents'""")
    
    parser.add_argument("--species_list_file", default="species_list.txt", help="File name to save the list of species names when --list_all_species is enabled.")
    parser.add_argument("-d", "--download_species", action='store_true', default=False,
                        help="Download species available in NCBI, limited by --search_term or --accepted_divisions and ")
    parser.add_argument("-m", "--max_species", type=int, default=10000,
                        help="Maximum number of species to download (default: 10k)")
    parser.add_argument("-ft", "--file_types", nargs="+", default=['genome', 'gtf'],
                        help="""List of assets to download e.g. genome gtf""")
    parser.add_argument("-e", "--email", type=str, default="bioinfo@kaust.edu.sa",
                        help="Email address for NCBI Entrez (default: bioinfo@kaust.edu.sa)")
    parser.add_argument("-o", "--output", type=str, default="genome_data",
                        help="Output directory (default: genome_data)")
    parser.add_argument("-p", "--processes", type=int, default=1,
                        help="Number of parallel downloads (default: 1)")
    args = parser.parse_args()
    
    #list all species
    if args.list_all_species:
        species = fetch_ncbi_species(args.email, args.max_species, args.accepted_divisions, args.species_list_file)
        print(f"Total species found: {len(species)}")
        for s in species[:10]:
            print(s)
            sys.exit(0)

    if args.search_term:
        print('Species found:\n', '\n'.join(get_species_names(search_term=args.search_term, accepted_divisions=args.accepted_divisions)))
        if not args.download_species:
            sys.exit(0)
    
    if args.max_species>10 and (args.search_term or args.list_all_species):
        print(f"Retrieving the list of {args.max_species} may take very long time..... It is faster to provide an input file with a long list of species. \nProcess Started... ")
    # Create output directory
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Get list of species either from command line or file
    species_list = []
    if args.species:
        species_list = args.species
    elif args.file and os.path.isfile(args.file):
        species_list = [line.strip() for line in open(args.file, 'r').readlines() if line.strip() if line!=""]
    elif args.search_term and args.download_species:
        species_list = get_species_names(search_term=args.search_term, accepted_divisions=args.accepted_divisions)
    elif args.download_species and args.accepted_divisions!=['all']:
        species_list = fetch_ncbi_species(args.email, args.max_species, args.accepted_divisions, args.species_list_file)
    else:
        print("Please provide a list of species or a file containing species names.")
        sys.exit(1)

    # Process species in parallel
    successful = 0
    failed = 0
    failed_species = []
    
    results = []
    if args.processes>1:
        with Pool(processes=args.processes) as pool:
            results = pool.starmap(download_genome, [(spcies, output_dir, ','.join(args.file_types)) for spcies in species_list])
    else:
        for spcies in species_list:
            results.append(download_genome(spcies, output_dir, ','.join(args.file_types)))
    print(results)

if __name__ == "__main__":
    main()
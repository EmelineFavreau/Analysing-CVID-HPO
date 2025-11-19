# --------------------------------------------------------
# 1. Setup and Configuration
# --------------------------------------------------------
import pandas as pd
import networkx as nx
from node2vec import Node2Vec
import argparse 

# --------------------------------------------------------
# 2. Argument Parsing Function
# --------------------------------------------------------
def parse_args():
    """Parses command line arguments for input and output file paths."""
    parser = argparse.ArgumentParser(
        description="Generate Node2Vec embeddings from a tab-separated edgelist."
    )
    # Define the argument for the input TSV file
    parser.add_argument(
        'input_file', 
        type=str, 
        help="Path to the input TSV file containing the edgelist (Source, Target)."
    )
    # Define the argument for the output embedding file
    parser.add_argument(
        'output_file', 
        type=str, 
        help="Path for the output Node2Vec embedding file (.emb format)."
    )
    return parser.parse_args()

# --------------------------------------------------------
# 3. Main Execution Function
# --------------------------------------------------------
def main():
    # Parse the arguments
    args = parse_args()
    INPUT_FILENAME = args.input_file
    EMBEDDING_FILENAME = args.output_file

    # Load Data and Build Graph
    try:
        # Assuming the file is tab-separated and contains two columns (Source and Target)
        df = pd.read_table(INPUT_FILENAME, header=None)
    except FileNotFoundError:
        print(f"Error: Input file '{INPUT_FILENAME}' not found. Exiting.")
        return

    # Create an undirected graph from the edgelist.
    G = nx.from_pandas_edgelist(df, source=0, target=1, create_using=nx.Graph())

    # Generate Embeddings with Node2Vec
    node2vec = Node2Vec(
        G,
        dimensions=128,
        walk_length=10,
        num_walks=80,
        workers=4,
        p=1,
        q=1
    )

    # Embed nodes using gensim's Word2Vec model
    model = node2vec.fit(window=10, min_count=1, batch_words=4)

    # Save Embeddings
    model.wv.save_word2vec_format(EMBEDDING_FILENAME)

if __name__ == '__main__':
    main()
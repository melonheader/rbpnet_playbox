import pandas as pd

def process_fasta(filepath: str, inject='identity', storage='dict', output_file=None, match_name=None, match_seq=None, **inject_kwargs):
    """Process a FASTA file entry by entry, apply a function to the sequences, and store or write the results.

    Args:
        filepath (str): Path to the FASTA file.
        inject (str, optional): Function to apply to each sequence, or 'identity' to return the sequence unchanged. The function can return either a string or a dictionary. Defaults to 'identity'.
        storage (str, optional): 'dict' to store results in a dictionary, 'file' to write results to a FASTA file, or 'dataframe' to store results in a pandas DataFrame. Defaults to 'dict'.
        output_file (str, optional): Path to the output file if storage is set to 'file' or 'dataframe'. Defaults to None.
        match_name (str, optional): A string to match the entry names. Only entries with names containing this string are processed. Defaults to None.
        match_seq (str, optional): A string to match the sequences. Only entries with sequences containing this string are processed. Defaults to None.
        inject_kwargs (dict, optional): Additional keyword arguments to pass to the inject function.

    Raises:
        ValueError: If an invalid inject function is provided.
        ValueError: If an invalid storage option is chosen, or if no output file is specified when storage is set to 'file'.

    Returns:
        dict or pandas.DataFrame: The processed results, either as a dictionary or a pandas DataFrame, depending on the storage option.
    """
    # initialise inject function
    if inject == 'identity':
        def inject_function(x):
            return x
    elif callable(inject):
        inject_function = inject
    else:
        raise ValueError("Invalid inject function. Provide a callable object")
    # initialize storage
    if storage == 'dict':
        results = {}
    elif storage == 'dataframe':
        rows = []
    elif storage == 'file' and output_file:
        storage = 'file'
    elif storage == 'file' and output_file is None:
        raise ValueError("An output file must be specified when storage is set to 'file'.")
    else:
        raise ValueError("Invalid storage option. Choose 'dict', 'file', or 'dataframe'.")
    # open the FASTA file connection
    with open(filepath, 'r') as fasta_file:
        # initialise
        entry_name = None
        sequence = []
        for line in fasta_file:
            line = line.strip()
            if line.startswith(">"):
                if entry_name:
                    full_sequence = ''.join(sequence)
                    if match_seq is None or match_seq in full_sequence:
                        processed_result = inject_function(full_sequence, **inject_kwargs)
                        if match_name is None or match_name in entry_name:
                            if storage == 'dict':
                                results[entry_name] = processed_result
                            elif storage == 'file':
                                with open(output_file, 'a') as out_file:
                                    if isinstance(processed_result, dict):
                                        for key, value in processed_result.items():
                                            out_file.write(f">{entry_name}|{key}\n{value}\n")
                                    else:
                                        out_file.write(f">{entry_name}\n{processed_result}\n")
                            elif storage == 'dataframe':
                                if isinstance(processed_result, dict):
                                    row = {'Entry_Name': entry_name}
                                    row.update(processed_result)
                                    rows.append(row)
                                else:
                                    rows.append({'Entry_Name': entry_name, 'Sequence': processed_result})
                # start an entry
                entry_name = line[1:]  # remove the ">"
                sequence = []
            else:
                # accumulate
                sequence.append(line)
        # finish the last entry
        if entry_name:
            full_sequence = ''.join(sequence)
            if match_seq is None or match_seq in full_sequence:
                processed_result = inject_function(full_sequence, **inject_kwargs)
                if match_name is None or match_name in entry_name:
                    if storage == 'dict':
                        results[entry_name] = processed_result
                    elif storage == 'file':
                        with open(output_file, 'a') as out_file:
                            if isinstance(processed_result, dict):
                                for key, value in processed_result.items():
                                    out_file.write(f">{entry_name}|{key}\n{value}\n")
                            else:
                                out_file.write(f">{entry_name}\n{processed_result}\n")
                    elif storage == 'dataframe':
                        if isinstance(processed_result, dict):
                            row = {'Entry_Name': entry_name}
                            row.update(processed_result)
                            rows.append(row)
                        else:
                            rows.append({'Entry_Name': entry_name, 'Sequence': processed_result})
    # handle storage
    if storage == 'dict':
        return results
    elif storage == 'dataframe':
        df = pd.DataFrame(rows)
        if output_file:
            df.to_csv(output_file, index=False)
        else:
            return df
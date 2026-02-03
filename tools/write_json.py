# This is a utility that attempts to write a dictionary to a
# file that is a little nicer than "json.dump()"
# While that utility does nothing wrong, I just dont
# like how it writes out vector data ie anything inside []
# on separate lines.  This utility will write anything inside an
# innermost bracket on one single line.

# ALSO: This is recursive, so the calling function should have last_att = True

import code  # For development: code.interact(local=dict(globals(), **locals()))
import json

def traverse_data(outfile, data: dict, indent_level: int = 0, current_key: str = "", last_att: bool = True) -> None:
    """
    Recursively traverses a nested dictionary or list and prints its contents.

    Args:
        data: The dictionary, list, or primitive value to process.
        indent_level: The current depth of recursion, used for formatting output.
        current_key: The key associated with the data (for dictionaries).
    """
    # Create the indentation string
    indent = "  " * indent_level

    # --- Case 1: The data is a dictionary (like the root or nested objects) ---
    if isinstance(data, dict):
        # Print the current dictionary key if one exists
        if current_key:
            outfile.write(f'{indent}"{current_key}":{{\n')
        else:
            outfile.write(f"{{\n")
        
        # Loop through all key-value pairs in the dictionary
        n_items = len(data.items())
        for i, (key, value) in enumerate(data.items()):
            is_last = (i == n_items - 1)
            # Recursively call the function for the value
            traverse_data(outfile,value, indent_level + 1, key, last_att = is_last)

        if(last_att):
            outfile.write(f"{indent}}}\n")
        else:
            outfile.write(f"{indent}}},\n")
            
    # --- Case 2: The data is a list (like an array in JSON) ---
    elif isinstance(data, list):

        for i,item in enumerate(data):
            if item is None:
                data[i] = "null"
                
        data_str = json.dumps(data)
        cleaned_data_str =  data_str.replace('"null"', 'null')

        if(last_att):
            outfile.write(f'{indent}"{current_key}": '+cleaned_data_str+"\n")
        else:
            outfile.write(f'{indent}"{current_key}": '+cleaned_data_str+",\n")
        
    # --- Case 3: The data is a primitive value (string, number, boolean, None) ---
    else:
        # Format the output for a primitive value
        value_type = type(data).__name__
        # Use repr() for strings to show quotes, helping differentiate values
        #display_value = repr(data) if isinstance(data, str) else str(data)
        display_value = '"'+data.strip()+'"' if isinstance(data,str) else str(data)
        if(last_att):
            outfile.write(f'{indent}"{current_key}": {display_value}\n')
        else:
            outfile.write(f'{indent}"{current_key}": {display_value},\n')


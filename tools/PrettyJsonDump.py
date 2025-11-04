# This is a utility that attempts to write a dictionary to a
# file that is a little nicer than "json.dump()"
# While that utility does nothing wrong, I just dont
# like how it writes out vector data ie anything inside []
# on separate lines.  This utility will write anything inside an
# innermost bracket on one single line.

#def PrettyJSONDump(json_dict,filename):

def traverse_data(data: ComplexData, indent_level: int = 0, current_key: str = "") -> None:
    """
    Recursively traverses a nested dictionary or list and prints its contents.

    Args:
        data: The dictionary, list, or primitive value to process.
        indent_level: The current depth of recursion, used for formatting output.
        current_key: The key associated with the data (for dictionaries).
    """
    # Create the indentation string
    indent = "    " * indent_level

    # --- Case 1: The data is a dictionary (like the root or nested objects) ---
    if isinstance(data, dict):
        # Print the current dictionary key if one exists
        if current_key:
            print(f"{indent}**Dictionary Key:** {current_key} (Size: {len(data)})")
        else:
            print(f"{indent}**Starting Traversal of Root Dictionary** (Size: {len(data)})")
        
        # Loop through all key-value pairs in the dictionary
        for key, value in data.items():
            # Recursively call the function for the value
            traverse_data(value, indent_level + 1, key)

    # --- Case 2: The data is a list (like an array in JSON) ---
    elif isinstance(data, list):
        print(f"{indent}**List Key:** {current_key} (Items: {len(data)})")
        # Loop through all items in the list using an index
        for index, item in enumerate(data):
            # Recursively call the function for the item
            # We use the key name plus the index for better context
            traverse_data(item, indent_level + 1, f"{current_key}[{index}]")

    # --- Case 3: The data is a primitive value (string, number, boolean, None) ---
    else:
        # Format the output for a primitive value
        value_type = type(data).__name__
        # Use repr() for strings to show quotes, helping differentiate values
        display_value = repr(data) if isinstance(data, str) else str(data)
        
        print(f"{indent}**Value Key:** {current_key} | **Type:** {value_type} | **Value:** {display_value}")


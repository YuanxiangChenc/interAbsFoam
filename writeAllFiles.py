import os
from pathlib import Path

def combine_files_in_current_folder(output_filename="combined_output.txt"):
    # Target the current working directory
    current_dir = Path(".")
    output_path = current_dir / output_filename

    # Open the final output file for writing (using UTF-8 to prevent encoding issues)
    with open(output_path, "w", encoding="utf-8") as outfile:
        # Loop through all items in the directory
        for item in current_dir.iterdir():
            # Process only files, skip directories, and ignore the output file itself
            if item.is_file() and item.name != output_filename:
                try:
                    # Write a visual header indicating which file's content follows
                    outfile.write(f"\n--- START OF FILE: {item.name} ---\n")
                    
                    # Read the contents of the file and write to the combined file
                    with open(item, "r", encoding="utf-8", errors="ignore") as infile:
                        outfile.write(infile.read())
                        
                    outfile.write(f"\n--- END OF FILE: {item.name} ---\n")
                    print(f"Successfully merged: {item.name}")
                    
                except Exception as e:
                    print(f"Could not read {item.name} due to error: {e}")

    print(f"\nFinished! All contents saved to: {output_path.resolve()}")

if __name__ == "__main__":
    combine_files_in_current_folder()


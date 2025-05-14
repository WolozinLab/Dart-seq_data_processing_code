#!/usr/bin/env python3
import os
import sys
import subprocess
import datetime

def get_latest_commit_message():
    """Get the latest git commit message"""
    try:
        result = subprocess.run(
            ["git", "log", "-1", "--pretty=%B"], 
            capture_output=True, 
            text=True, 
            check=True
        )
        return result.stdout.strip()
    except subprocess.CalledProcessError:
        return "Unknown commit"

def update_readme():
    """Update the README.md file with the latest commit information"""
    readme_path = "README.md"
    
    # Check if README.md exists
    if not os.path.exists(readme_path):
        print(f"Error: {readme_path} not found")
        return 1
    
    # Read the current README content
    with open(readme_path, 'r') as file:
        content = file.readlines()
    
    # Find the Changelog section
    changelog_index = -1
    for i, line in enumerate(content):
        if line.strip() == "## Changelog":
            changelog_index = i
            break
    
    if changelog_index == -1:
        print("Error: Changelog section not found in README.md")
        return 1
    
    # Get current date and commit message
    current_date = datetime.datetime.now().strftime("%Y-%m-%d")
    commit_message = get_latest_commit_message()
    
    # Create new changelog entry
    new_entry = f"- {current_date}: {commit_message}\n"
    
    # Insert new entry after the Changelog heading
    content.insert(changelog_index + 1, new_entry)
    
    # Write updated content back to README.md
    with open(readme_path, 'w') as file:
        file.writelines(content)
    
    print(f"Updated README.md with new changelog entry: {new_entry.strip()}")
    return 0

if __name__ == "__main__":
    sys.exit(update_readme()) 
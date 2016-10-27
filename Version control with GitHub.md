# Purpose
Putting your existing work on GitHub and make version control

# Step-by-step pipeline for a new repository
1. Create a new repository on GitHub. To avoid errors, do not initialize the new repository with README, license, or gitignore files. 

2. Open Terminal and change the current working directory to your local project.

3. Initialize the local directory as a Git repository.
`$ git init`

4. Add the files in your new local repository and stage them for commit.
`$ git add .`

5. Commit the files that you’ve staged in your local repository.
`$ git commit -m “First commit”`

6. At the top of your GitHub repository’s Quick Setup page, click  to copy the remote repository URL.

7. Set the `origin` with the URL you got in step 6 with command `git remote add origin https://github.com/......`

7. In Terminal, add the URL for the remote repository where your local repository will be pushed.
`$ git push origin master`

# Make change or adding a file to an existing repository 
1. On your computer, make change in an existing file or move a new file to into the local directory.

2. Open Terminal and change the current working directory to your local project.

3. Stage the file for the first commit to your repository.
`$ git add .`

4. Commit the files that you’ve staged in your local repository.
`$ git commit -m “Adding commit”`

5. Push the changes in your local repository to GitHub.
`$ git push origin master`

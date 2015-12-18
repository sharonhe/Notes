# Basic Unix Commands


- Open a new tab in Terminal: `Command + T`

- Open Markdown in browser: `Command + shift + p`

- :heart_eyes: Stop running: `Ctrl + Z`

- Exit to command prompt: `Ctrl + C`

- Go to parent directory: `cd ..`

- Go to directory on a external hard drive, eg AdyBox
`cd Volumnes/AdyBox/`

- Print work directory
`pwd`

- List files with detailed information
`ls -l`

- List files with detailed information, sort by modification time
`ls -lt`

- List files by incomplete file name
```
ls apple.*
ls ?pple.genome
ls [a-h]*.genome
ls {pear,peach}.genome
```

- Look at mannual pages
`man ls`

- Make directories
`mkdir apple pear peach`

- Copy files to a directory
`cp apple.genome apple.sample apple`

- Move files to a directory
`mv apple.genome apple.sample apple`

- Remove a file from current directory
`rm apple.genome`

- Tend to Remove a file from current directory, a questioned confirmation will be needed before excution
`rm -i apple.genome`

- Remove empty directory
`rmdir apple`

- Remove all the files in a directory and the directory itself
`rm -r apple/`

- Access content of a file
`more apple.genome`
`less apple.genome`
`head -50 apple.genome`
`tail -15 apple.genome`

- Look at the next chromosome sequence in a FASTA file
`/>`

- Concatenate multiple files
`cat peach.genes apple.genes`
`cat */*.genes`

- Concatenate multiple files, and allow us to see a page at a time
`cat */*.genome | more`

- Word count
`wc apple.genome`
`wc -l apple.genome`

- Redirect an output from terminal into a file
`wc -l apple.genome > nlines.wc`

- Redirect a result from an input file to the terminal
`wc -l < apple.genome`

- Make the output to be an input of another command
`ls apple.genes | wc -l`

- Sort files
`sort months`
`sort -r months`
`sort -k 2 months`
`sort -k 2n months`
`sort -k 2nr months`
`sort -k 3 -k 2n months`

- Extract Columne(s) separated by tabular
`cut -f1 months`
`cut -f1,2 months`
`cut -f1-3 months`

- Extract Columne(s) separated by space
`cut -d ' ' -f1 months`

- Sort uniquely
`sort -u seasons`
`sort seasons | uniq`

- Sort uniquely, and show how many repeats for each word
`sort seasons | uniq -c`

- Find the name content in multiple files
`grep "root" */*.samples`

- Find the name content in multiple files, and show the line number in the file contains this content
`grep -n "root" */*.samples`

- Compress/uncompress one file
`gzip apple.genome`
`gunzip apple.genome.gz`
`bzip2 apple.genome`
`bunzip2 apple.genome.bz2`

- Archive multiple files, and compress it
`tar -cvf Apple.tar apple.genes apple.genome apple.samples`
`gzip Apple.tar`

- Uncompress the file and Unarchive multiple files
`gunzip Apple.tar.gz`
`tar -xvf Apple.tar`

# Sublime
## Column selection across the entire file
1. Ctrl+A - select all.
2. Ctrl+Shift+L - split selection into lines.
3. Then move all cursors with left/right, select with Shift+left/right. Move all cursors to start of line with Home.



























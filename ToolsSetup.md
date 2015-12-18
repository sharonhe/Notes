
# How to setup tools on OS X

## How to install XCode Command Line Tools

- `xcode-select --install`

## How to install `wget`

- Install HomeBrew form here: http://brew.sh/

- `ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"`

- `brew install wget`

## How to install and configure SRA Toolkit

- SRA Toolkit Installation and Configuration Guide can be found here: http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc&f=std

## How to set `$PATH`

- `echo $PATH`

- Should look like `/usr/bin:/bin:/usr/sbin:/sbin:/usr/local/bin`

- To add a permanent location to `$PATH`, follow these steps

  - `cd`
  - `nano .bash_profile`
  - This will open nano editor. Add this line to the file `export PATH=“/your/path1:/your/path2:/your/path3:$PATH"`

  ![Image](http://coolestguidesontheplanet.com/wp-content/uploads/2014/02/osx-path-modify1.png)


  - ‘control’ +’o’ and Return to Save
  - ‘control’ +’x’ to exit nano
  - Restart terminal
  - `echo $PATH` 


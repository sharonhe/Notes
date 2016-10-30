# How to install IGVtools

1. Install IGVTools from https://software.broadinstitute.org/software/igv/download (Not the GUI, just the commandline tools)
2. CD to that folder and try to run it from command line with `./igvtools`
3. If it complains about java major.minor version, you probably have the wrong java version. You can check it with `java -version`
4. You need to update to JDK 1.7 or 1.8. Download it from [here](http://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html).
5. Install it by double clicking the `.dmg` file
6. Now go back to terminal and check `java -version`. If it still points to the old version, you can point your system to the latest version you installed with 
```
export JAVA_HOME=`/usr/libexec/java_home -v 1.8`
```
And then recheck your `java -version`
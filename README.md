# Raj Lab single molecule RNA FISH probe design software
Largely written in MATLAB, sorry.


As of August 2024, connection with Repeat Masker does not work. 
In order to use Repeat Masker manually, you can upload your fasta file to RepeatMasker's website and download the corresponding masked file.
Then, you can use findprobesLocal, setting 'repeatmask' and 'repeatmaskmanual' to True and setting 'repeatmaskfile' to the path to your file.
For example:
`findprobesLocal('mCherry.fasta',32,'repeatmask', true,'repeatmaskmanual',true,'repeatmaskfile','mCherry_repeatmasked.fa','species','mouse')`

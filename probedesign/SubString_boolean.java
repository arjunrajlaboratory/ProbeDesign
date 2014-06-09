import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

public class SubString_boolean {
	private HashSet<String> geneSet = new HashSet<String>(); //contains substrings of a gene
	private boolean[] found; //has substring been found in database?
	private StringBuilder outseq;
	private long findtic; private long findtoc; private long totalSearchTime=0;

    // ** CONSTRUCTOR METHOD ** //
	public SubString_boolean(String filename, int subLength, String inseq) throws IOException {
		/* Input
		 * filename - full path of the database FASTA file to be searched (can be multi sequence)
		 * subLength - length of substring queries
		 * inseq - query sequence (eg. the mRNA sequence that you want to mask against the database)
		 * */
    
        /* Usage in MATLAB
        * >> javaaddpath('/Users/marshall/rajlab/sequenceanalysis/probedesign/')
        * >> genereader = SubString_boolean('/path/to/probedesign/pseudogeneDBs/human.fasta', 16, upper(inseq));
        * >> maskedSeq = char(genereader.FindSeq());  % convert java.String to MATLAB char()
        */
		
		// boolean array with same length as query sequence minus sub-string query length
		found = new boolean[inseq.length()-(subLength-1)];  // default is all FALSE
         
        // make a copy of "inseq" that will be edited later with "X" for match positions
        outseq = new StringBuilder(inseq);

        // load the database file into a BufferedReader
		if (filename == null) { throw new IllegalArgumentException(); }
		BufferedReader reader = new BufferedReader (new FileReader(filename));

		/* Logic for reading and searching database file: 
         *  1. We read the file line-by-line until we complete a single sequence in a (potentially)
         *     multi-sequence FASTA file. 
         *  2. Break up the single gene sequence into substrings of length=subLength and add
         *     them to the HashSet geneSet. This is a unique collection of substrings.
         *  3. Perform search for substrings of inseq in HashSet geneSet. Positions of matching 
         *     substrings in inseq are recorded in the boolean array "found" as TRUE
         *  4. Clear the HashSet "geneSet" and repeat steps 1-3.   
         */
		if (reader.ready()){
			String newLine = reader.readLine().trim(); //first line to read; runs only once
			String newGene = new String(); //initialize container for one gene

			while(newLine != null){ //while there are still lines to read

				if (newLine.contains(">")){ //if new line is a gene title
					newGene = new String(); //clear gene container
					newLine = reader.readLine().trim(); //read next line
					continue; //skip code below this line
				} else {

					while (!newLine.contains(">")){ //while new line is a continuation of the gene
						newLine = newLine.trim();
						newGene = newGene.concat(newLine); //add new line to gene container
						newLine = reader.readLine(); //read next line
						if (newLine == null){break;} //stop at the end of the database
					}

					//enters substring entries of a single FASTA sequence into HashSet geneSet
					for (int i=0; i<newGene.length()-(subLength-1); i++){ 
						String subString = newGene.substring(i,i+subLength); //generates hashcode from substring
						geneSet.add(subString);
					}
					
                    // perform search for substrings of inseq in HashSet geneSet
					findtic = System.currentTimeMillis();
					for (int i=0; i<inseq.length()-(subLength-1); i++){ //for length of query sequence
						if (!found[i]){ //if already found, search can be skipped
							found[i] = geneSet.contains(inseq.substring(i,i+subLength)); //was sequence found in this segment of genes?
							}
						}
					findtoc = System.currentTimeMillis();
					totalSearchTime += findtoc-findtic; //pools HashSet search time
					
					geneSet.clear(); //clear HashSet and process next gene
					
				}
			}	
		}
		reader.close();
		
		System.out.println("Total HashSet search time: " + Long.toString(totalSearchTime/1000) + " seconds"); //outputs total search time
		
		// Make a copy of "inseq" and replace characters that had matches to database with "X"
		for (int i=0; i<inseq.length()-(subLength-1); i++){
			if (found[i]){
				outseq.replace(i,i+1,"X");
			}
		}
	}

    // ** METHOD TO RETURN OUTPUT: your input sequence masked with X's ** //
	public String FindSeq(){
		return outseq.toString();
	}
}

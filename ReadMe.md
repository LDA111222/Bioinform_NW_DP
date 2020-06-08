Student Name: WangQi 王琦.
Student Id: 1710476.

This program uses Smith-Waterman algorithms.
The two sequences should save as "*.txt".(ex. Test1.txt or Test2.txt)
All the letters are excepted to be closely connected without spaces, line breaks or other extraneous characters.
Uppercase, lowercase and mixture are allowed.(ex. AAACT aaact aAAcT)
The score matrix should save as "*.txt".(ex. BLOSUM62.txt)
Warning: All letters in the sequence files should appear in the score matrix.
The format should be consistent with the file BLOSUM62.txt(Download from http://yanglab.nankai.edu.cn/teaching/bioinformatics/BLOSUM62.txt)
JDK(version >= 8) is necessary, although older version may work.(Test Program runs with java version "1.8.0_202", Java(TM) SE Runtime Environment (build 1.8.0_202-b08), Java HotSpot(TM) 64-Bit Server VM (build 25.202-b08, mixed mode))
To run the program,
1. cd in the dirname
2. javac Practice2.java
3. java Practice2 ScoreMatrixFile.txt Sequence1File.txt Sequence2File.txt
The output result is composed of three part.
The default value of a and b are -3 and -1, to modify a, b please open the Practice2.java, in Line 390(int a = -3, b = -1;), key in the value after '='.
The first part is the input sequence comfirmation.
The second part is the matrix M, matrix Iy and matrix Ix.(Numbers which are no more than -10000000 means negative infinity.)
The last part is the alignment result.

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;

class ScoreMatrix {
    protected int N;
    protected ArrayList<Integer> dataOfMatrix = new ArrayList<Integer>();
    protected ArrayList<Character> aminoSequence = new ArrayList<Character>();

    void readMatrix(String filename) {
        try {
            FileReader input = new FileReader(filename);
            BufferedReader bufInput = new BufferedReader(input);
            int ch, sentry = 0, sym = 1, buck = 0;
            //
            while ((ch = bufInput.read()) != -1) {

                if ((ch >= 65 && ch <= 90 || ch == 42) && buck == 0) {
                    aminoSequence.add((char) ch);

                    if (ch == 42) {
                        buck = 1;
                        this.N = aminoSequence.size();
                    }

                }

                if (ch >= 48 && ch <= 57 && sentry == 0) {
                    dataOfMatrix.add(sym * (ch - 48));
                    sym = 1;
                    sentry = 1;
                } else if (sentry == 1 && ch >= 48 && ch <= 57) {
                    dataOfMatrix.set(dataOfMatrix.size() - 1,
                            dataOfMatrix.get(dataOfMatrix.size() - 1) * 10 + (ch - 48));
                    sym = 1;
                    sentry = 0;
                } else if (ch == '-') {
                    sym = -1;
                } else {
                    sentry = 0;
                }
            }

            input.close();
            bufInput.close();
        } catch (IOException e1) {
            e1.printStackTrace();
        }
    }

    void printMatrix() {

        for (int j = 0; j < aminoSequence.size(); j++) {
            System.out.print(aminoSequence.get(j) + " ");
        }

        for (int j = 0; j < dataOfMatrix.size(); j++) {
            System.out.print(dataOfMatrix.get(j) + " ");
            if ((j + 1) % 24 == 0) {
                System.out.print('\n');
            }
        }

        System.out.printf("The dim of Matrix is %d^2\nThe size of Matrix is %d", aminoSequence.size(),
                dataOfMatrix.size());

    }

    ScoreMatrix(String filename) {
        readMatrix(filename);
    }
}

class analySequence {
    protected int N;
    protected ArrayList<Character> analysNeed = new ArrayList<Character>();

    analySequence(String filename) {
        try {
            FileReader input = new FileReader(filename);
            BufferedReader bufInput = new BufferedReader(input);
            int ch;

            while ((ch = bufInput.read()) != -1) {

                if (ch >= 65 && ch <= 90 || ch == 42) {
                    analysNeed.add((char) ch);
                } else if (ch >= 97 && ch <= 122) {
                    analysNeed.add((char) (ch - 32));
                }

            }
            input.close();
            bufInput.close();
            N = analysNeed.size();
        } catch (IOException e1) {
            e1.printStackTrace();
        }
    }

    void printAnalySequence() {

        for (int j = 0; j < analysNeed.size(); j++) {
            System.out.print(analysNeed.get(j) + " ");
        }
        System.out.printf("\nThe sequence length is %d", N);
    }
}

class globalCompare {
    protected int MAX_M = 1000, MAX_N = 1000;
    protected int N = 100, M = 100, a, b, s1=1, s2=-1;
    protected ArrayList<Integer> dataOfMatrix = new ArrayList<Integer>();
    protected ArrayList<Character> aminoSequence = new ArrayList<Character>();
    // N is the length of seq1
    // M is the length of seq2
    // a and b are the fine score when it is extended
    // s1 is the M scores when the animo acids were the same while s2 is the score
    // when they're not
    protected ArrayList<Character> seq1 = new ArrayList<Character>();
    protected ArrayList<Character> seq2 = new ArrayList<Character>();
    protected int[][] processMatrixM = new int[MAX_N + 1][MAX_M + 1], processMatrixIx = new int[MAX_N + 1][MAX_M + 1],
            processMatrixIy = new int[MAX_N + 1][MAX_M + 1];
    protected int[][] pathMatrixM = new int[MAX_N + 1][MAX_M + 1], pathMatrixIx = new int[MAX_N + 1][MAX_M + 1],
            pathMatrixIy = new int[MAX_N + 1][MAX_M + 1];
    // The number in the pathMatrix: 0 means from i-1, j-1; -1 means from i, j-1; 1 means from i-1, j. 
    globalCompare(ArrayList<Character> l1, ArrayList<Character> l2, int a, int b, int N, int M, ArrayList<Integer> dataOfMatrix, ArrayList<Character> aminoSequence) {
        this.seq1 = l1;
        this.seq2 = l2; 
        this.dataOfMatrix = dataOfMatrix;
        this.aminoSequence = aminoSequence;
        this.a = a;
        this.b = b;
        this.N = N;
        this.M = M;
        this.processMatrixM[0][0] = 0;
        this.processMatrixIx[0][0] = a;
        this.processMatrixIy[0][0] = a;
        this.pathMatrixM[0][0] = 2;
        this.pathMatrixIx[0][0] = 2;
        this.pathMatrixIy[0][0] = 2;

        for (int i = 1; i <= M; i++) {
            this.processMatrixM[0][i] = -10000000; // -10e7
            this.processMatrixIx[0][i] = -10000000; // -10e7
            this.processMatrixIy[0][i] = a + b * (i);
            this.pathMatrixIy[0][i] = -1;
        }

        for (int i = 1; i <= N; i++) {
            this.processMatrixM[i][0] = -10000000; // -10e7
            this.processMatrixIx[i][0] = a + b * (i);
            this.processMatrixIy[i][0] = -10000000; // -10e7
            this.pathMatrixIx[i][0] = 1;
        }
        // printamino();
    }
    globalCompare(int a, int b, ArrayList<Integer> dataOfMatrix, ArrayList<Character> aminoSequence) {
        this.dataOfMatrix = dataOfMatrix;
        this.aminoSequence = aminoSequence;
        this.a = a;
        this.b = b;
        this.N = 100;
        this.M = 100;
        this.processMatrixM[0][0] = 0;
        this.processMatrixIx[0][0] = a;
        this.processMatrixIy[0][0] = a;
        this.pathMatrixM[0][0] = 2;
        this.pathMatrixIx[0][0] = 2;
        this.pathMatrixIy[0][0] = 2;
        randomMatrixGet();

        for (int i = 1; i <= M; i++) {
            this.processMatrixM[0][i] = -10000000; // -10e7
            this.processMatrixIx[0][i] = -10000000; // -10e7
            this.processMatrixIy[0][i] = a + b * (i);
            this.pathMatrixIy[0][i] = -1;
        }

        for (int i = 1; i <= N; i++) {
            this.processMatrixM[i][0] = -10000000; // -10e7
            this.processMatrixIx[i][0] = a + b * (i);
            this.processMatrixIy[i][0] = -10000000; // -10e7
            this.pathMatrixIx[i][0] = 1;
        }
    }
    void randomMatrixGet(){
        System.out.println(aminoSequence.size());
        for (int i = 0; i < 100; i++){
                int rand1 = (int)Math.floor(Math.random()*23), rand2 = (int)Math.floor(Math.random()*23);
                seq1.add(aminoSequence.get(rand1));
                seq2.add(aminoSequence.get(rand2));
        }
        for (int i = 0; i < seq1.size(); i++) {
            System.out.print(seq1.get(i));
        }
        System.out.println();
        System.out.println(seq1.size());
        for (int i = 0; i < seq2.size(); i++) {
            System.out.print(seq2.get(i));
        }
        System.out.println();
        System.out.println(seq2.size());
    }
    /* void printamino(){
        System.out.println(aminoSequence.size());
        for (int i =0; i < aminoSequence.size(); i++ ){
            System.out.println(aminoSequence.get(i));
        }
    } */
    int ScoreSeek(char a, char b) {
        int a_num = aminoSequence.indexOf(a), b_num = aminoSequence.indexOf(b);
        if(a_num != -1 && b_num != -1) {
            return dataOfMatrix.get(a_num*aminoSequence.size()+b_num);
        }
        else {
            return -1;
        }
    }        

    void mState(int i, int j) {
        int sOfState = ScoreSeek(seq1.get(i-1), seq2.get(j-1));

        // if (seq1.get(i-1) == seq2.get(j-1)) {
        //     sOfState = s1;
        // } else {
        //     sOfState = s2;
        // }
        

        if (processMatrixM[i - 1][j - 1] > processMatrixIx[i - 1][j - 1]
                && processMatrixM[i - 1][j - 1] > processMatrixIy[i - 1][j - 1]) {
            this.pathMatrixM[i][j] = 0;
            this.processMatrixM[i][j] = (sOfState+processMatrixM[i - 1][j - 1]);
        } else if (processMatrixIx[i - 1][j - 1] > processMatrixIy[i - 1][j - 1]) {
            this.pathMatrixM[i][j] = 1;
            this.processMatrixM[i][j] = (sOfState+processMatrixIx[i - 1][j - 1]);
        } else {
            this.pathMatrixM[i][j] = -1;
            this.processMatrixM[i][j] = (sOfState+processMatrixIy[i - 1][j - 1]);
        }

    }

    void xState(int i, int j) {
        int sOfStateM = processMatrixM[i - 1][j] + a;
        int sOfStateIx = processMatrixIx[i - 1][j] + b;
        if (sOfStateM > sOfStateIx) {
            this.pathMatrixIx[i][j] = 0;
            this.processMatrixIx[i][j] = sOfStateM;
        } else {
            this.pathMatrixIx[i][j] = 1;
            this.processMatrixIx[i][j] = sOfStateIx;
        }
    }

    void yState(int i, int j) {
        int sOfStateM = processMatrixM[i][j - 1] + a;
        int sOfStateIy = processMatrixIy[i][j - 1] + b;
        if (sOfStateM > sOfStateIy) {
            this.pathMatrixIy[i][j] = 0;
            this.processMatrixIy[i][j] = sOfStateM;
        } else {
            this.pathMatrixIy[i][j] = -1;
            this.processMatrixIy[i][j] = sOfStateIy;
        }
    }

    void mainCompare() {

        for (int i = 1; i <= N; i++) {
            
            for (int j = 1; j <= M; j++) {
                mState(i, j);
                xState(i, j);
                yState(i, j);
            }

        }
        /* System.out.println("MatrixM");
        for (int j = 0; j <= M; j++) {
            
            for (int i = 0; i <= N; i++) {
               if(processMatrixM[i][j]<-100000){
                   System.out.print("-∞" + " ");
               }else{
                System.out.print(processMatrixM[i][j]+" ");
               }
            }
            System.out.println();

        }
        System.out.println("MatrixIx");
        for (int j = 0; j <= M; j++) {
            
            for (int i = 0; i <= N; i++) {
                if(processMatrixIx[i][j]<-100000){
                    System.out.print("-∞" + " ");
                }else{
                    System.out.print(processMatrixIx[i][j]+" ");
                }
            }
            System.out.println();

        }
        System.out.println("MatrixIy");
        for (int j = 0; j <= M; j++) {
            
            for (int i = 0; i <= N; i++) {
                if(processMatrixIy[i][j]<-100000){
                    System.out.print("-∞" + " ");
                }else{
                    System.out.print(processMatrixIy[i][j]+" ");
                }
            }
            System.out.println();

        }
        System.out.println("pathMatrixM");
        for (int j = 0; j <= M; j++) {
            
            for (int i = 0; i <= N; i++) {
               System.out.print(pathMatrixM[i][j]+" ");
            }
            System.out.println();

        }
        System.out.println("pathMatrixIx");
        for (int j = 0; j <= M; j++) {
            
            for (int i = 0; i <= N; i++) {
               System.out.print(pathMatrixIx[i][j]+" ");
            }
            System.out.println();

        }
        System.out.println("pathMatrixIy");
        for (int j = 0; j <= M; j++) {
            
            for (int i = 0; i <= N; i++) {
               System.out.print(pathMatrixIy[i][j]+" ");
            }
            System.out.println();

        }
        System.out.println(); */
    }
    void printResult() {
        int i = N, j = M, parot = 0;
        ArrayList<Character> resIx = new ArrayList<Character>();
        ArrayList<Character> resIy = new ArrayList<Character>();  
        if (this.processMatrixM[i][j] > this.processMatrixIx[i][j] && this.processMatrixM[i][j] > this.processMatrixIy[i][j]) {
            resIx.add(seq1.get(i-1));
            resIy.add(seq2.get(j-1));
            parot = pathMatrixM[i][j];
            i--;
            j--;
        }else if (this.processMatrixIx[i][j] >  this.processMatrixIy[i][j]) {
            resIx.add(seq1.get(i-1));
            resIy.add('-');
            parot = pathMatrixIx[i][j];
            i--;
        }else {
            resIx.add('-');
            resIy.add(seq2.get(j-1));
            parot = pathMatrixIy[i][j];
            j--;
        }
        while(i > 0 || j > 0){
            switch(parot){
                case 0:{
                    resIx.add(seq1.get(i-1));
                    resIy.add(seq2.get(j-1));
                    parot = pathMatrixM[i][j];
                    i--;
                    j--;
                    break;
                }
                case 1:{
                    resIx.add(seq1.get(i-1));
                    resIy.add('-');
                    parot = pathMatrixIx[i][j];
                    i--;
                    break;
                }
                case -1:{
                    resIx.add('-');
                    resIy.add(seq2.get(j-1));
                    parot = pathMatrixIy[i][j];
                    j--;
                    break;
                }
            }
        }
        /* while(i>1 && j>1) {    
            if (this.processMatrixM[i][j] > this.processMatrixIx[i][j] && this.processMatrixM[i][j] > this.processMatrixIy[i][j]) {
                resIx.add(seq1.get(i-1));
                resIy.add(seq2.get(j-1));
                switch (pathMatrixM[i][j]) {
                    case 0: {
                        i = i-1;
                        j = j-1;
                        // resIx.add(seq1.get(i-1));
                        // resIy.add(seq2.get(j-1));
                        break;
                    }
                    case 1: {
                        i = i-1;
                        // resIx.add(seq1.get(i-1));
                        // resIy.add('-');
                        break;
                    }
                    case -1: {
                        j = j-1;
                        // resIx.add('-');
                        // resIy.add(seq2.get(j-1));
                        break;
                    }
                }
            }else if (this.processMatrixIx[i][j] >  this.processMatrixIy[i][j]) {
                resIx.add(seq1.get(i-1));
                resIy.add('-');
                System.out.println(pathMatrixM[i][j]);
                switch (pathMatrixIx[i][j]) {
                    case 0: {
                        i = i-1;
                        j = j-1;
                        // resIx.add(seq1.get(i-1));
                        // resIy.add(seq2.get(j-1));
                        break;
                    }
                    case 1: {
                        j = j-1;
                        // resIx.add(seq1.get(i-1));
                        // resIy.add('-');
                        break;
                    }
                    case -1: {
                        i = i-1;
                        // resIx.add('-');
                        // resIy.add(seq2.get(j-1));
                        break;
                    }
                }
            }else {
                resIx.add('-');
                resIy.add(seq2.get(j-1));
                System.out.println(pathMatrixM[i][j]);
                switch (pathMatrixIy[i][j]) {
                    case 0: {
                        i = i-1;
                        j = j-1;
                        // resIx.add(seq1.get(i-1));
                        // resIy.add(seq2.get(j-1));
                        break;
                    }
                    case 1: {
                        i = i-1;
                        // resIx.add(seq1.get(i-1));
                        // resIy.add('-');
                        break;
                    }
                    case -1: {
                        j = j-1;
                        // resIx.add('-');
                        // resIy.add(seq2.get(j-1));
                        break;
                    }
                }
            }
        }
        // System.out.println("The length of res:" + resIx.size()); */
        System.out.println("I="+i+"J="+j);
        // while (i > 1) {
        //     resIx.add(seq1.get(i-1));
        //     resIy.add('-');
        //     i--;
        // }
        // while (j > 1) {
        //     resIx.add('-');
        //     resIy.add(seq2.get(j-1));
        //     j--;
        // }

        for (int ii = resIx.size()-1; ii >= 0; ii--) {
            System.out.print(resIx.get(ii));
        }
        System.out.println();
        for (int ii = resIy.size()-1; ii >= 0; ii--) {
            System.out.print(resIy.get(ii));
        }
    }


}

public class Practice2 {
    public static void main(String[] args) {
        int a = -11, b = -1;
        ScoreMatrix s1 = new ScoreMatrix(args[0]);
        analySequence a1 = new analySequence(args[1]);
        analySequence a2 = new analySequence(args[2]);
        System.out.println('\n');
        a1.printAnalySequence();
        System.out.println('\n');
        a2.printAnalySequence();
        System.out.println('\n');
        ArrayList<Character> l1 = a1.analysNeed;
        ArrayList<Character> l2 = a2.analysNeed;
        // globalCompare g1 = new globalCompare(a1.analysNeed, a2.analysNeed, -3, -1, a1.analysNeed.size(), a2.analysNeed.size(), s1.dataOfMatrix, s1.aminoSequence);
        globalCompare g1 = new globalCompare(l1, l2, a, b, l1.size(), l2.size(), s1.dataOfMatrix, s1.aminoSequence);
        // globalCompare g1 = new globalCompare(a, b, s1.dataOfMatrix, s1.aminoSequence);
        g1.mainCompare();
        g1.printResult();
        System.out.println();

    }
}

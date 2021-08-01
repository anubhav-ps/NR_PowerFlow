public class PowerFlow {

    //-------------variables---------------------
    // Magnitude of Voltage |Y|
    protected static float[][] Y_Matrix;
    // Theta_Matrix
    protected static float[][] Theta_Matrix;
    // bus complex load voltage
    protected static float[] Vol = new float[]{1, 1, 1.01f};  // Vol[0] -> V1
    // bus complex load angle
    protected static float[] Angle = new float[]{0, 0, 0};  // Angle[0] -> D1
    //active power demand of the 3 bus
    protected static float[] power = new float[]{0, -0.9f, 0.6f};
    //reactive power demand of the 3 bus
    protected static float[] q = new float[]{0, -0.5f, 0};
    //
    protected static float[] f = new float[]{0f, 0f, 0f};

    //-------------methods---------------------------
    //Constructor
    PowerFlow(float y1, float y2, float y3) {
        System.out.println("<<<<<<<<<<<<<<<<<<<< Power Flow Equation >>>>>>>>>>>>>>>>>>>>>>>>\n");
        System.out.println("\t\t\t\t\t\tDeveloped by : ");
        System.out.println("-----------------------------------------------------------------");
        System.out.printf("%-40s%-40s\n","Name","Register Number");
        System.out.println("-----------------------------------------------------------------");
        System.out.printf("%-40s%-40s\n","Aditya Vishwakarma","19BEE1178");
        System.out.printf("%-40s%-40s\n","Anubhav.P.S","19BEE1144");
        System.out.printf("%-40s%-40s\n","Chirag Rathore","19BEE1025");
        System.out.println("-----------------------------------------------------------------\n");
        System.out.println("<<<<<<<<<<<<<<<<<<<< Power Flow Equation >>>>>>>>>>>>>>>>>>>>>>>>\n");
        form_Y_Matrix(y1, y2, y3);
        form_Theta_Matrix();
        System.out.println("Entered Admittance value : ");
        System.out.printf("%-10s%-20f\n","Y1",y1);
        System.out.printf("%-10s%-20f\n","Y2",y2);
        System.out.printf("%-10s%-20f\n","Y3",y3);
    }

    // get the power flow equation
    public static float fOf(String eq) {
        double result = 0.0;
        if (eq.toCharArray()[0] == 'P') {  //if eq is active power demand
            if (eq.toCharArray()[1] == '1') {
                for (int i = 0; i < 3; i++) {
                    result += Y_Matrix[0][i] * Vol[i] * Vol[0] * HrMath.Cos(Theta_Matrix[0][i] + Angle[i] - Angle[0]);
                }
                result -= power[0];
                return HrMath.roundOff(result);
            } else if (eq.toCharArray()[1] == '2') {
                for (int i = 0; i < 3; i++) {
                    result += Y_Matrix[1][i] * Vol[i] * Vol[1] * HrMath.Cos(Theta_Matrix[1][i] + Angle[i] - Angle[1]);
                }
                result -= power[1];
                return HrMath.roundOff(result);
            } else if (eq.toCharArray()[1] == '3') {
                for (int i = 0; i < 3; i++) {
                    result += Y_Matrix[2][i] * Vol[i] * Vol[2] * HrMath.Cos(Theta_Matrix[2][i] + Angle[i] - Angle[2]);
                }
                result -= power[2];
                return HrMath.roundOff(result);
            }
        } else if (eq.toCharArray()[0] == 'Q') {  // if eq is reactive power demand
            if (eq.toCharArray()[1] == '1') {
                for (int i = 0; i < 3; i++) {
                    result += Y_Matrix[0][i] * Vol[i] * Vol[0] * HrMath.Sin(Theta_Matrix[0][i] + Angle[i] - Angle[0]);
                    result = result * (-1);
                }
                result -= q[0];
                return HrMath.roundOff(result);
            } else if (eq.toCharArray()[1] == '2') {
                for (int i = 0; i < 3; i++) {
                    result += Y_Matrix[1][i] * Vol[i] * Vol[1] * HrMath.Sin(Theta_Matrix[1][i] + Angle[i] - Angle[1]);
                    result = result * (-1);
                }
                result -= q[1];
                return HrMath.roundOff(result);
            } else if (eq.toCharArray()[1] == '3') {
                for (int i = 0; i < 3; i++) {
                    result += Y_Matrix[2][i] * Vol[i] * Vol[2] * HrMath.Sin(Theta_Matrix[2][i] + Angle[i] - Angle[2]);
                    result = result * (-1);
                }
                result -= q[2];
                return HrMath.roundOff(result);
            }
        }
        System.out.println("Invalid Input");
        return HrMath.roundOff(result);
    }

    //Y_Matrix
    private void form_Y_Matrix(float y1, float y2, float y3) {
        Y_Matrix = new float[][]{{y1 + y2, y1, y2}
                , {y1, y1 + y3, y3},
                {y2, y3, y3 + y2}};
    }

    private void form_Theta_Matrix() {
        float PI = (float) Math.PI;
        Theta_Matrix = new float[][]{{-PI / 2, PI / 2, PI / 2},
                {PI / 2, -PI / 2, PI / 2},
                {PI / 2, PI / 2, -PI / 2}};
    }

    public void display_Y_Matrix() {
        System.out.println("------------------Displaying Modulus of Y Matrix------------------");
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                float c1 = Y_Matrix[i][j];
                System.out.printf("%-28f", c1);
            }
            System.out.println();
        }
        System.out.println("---------------------------------------------------------------------\n\n");
    }

    public void display_Theta_Matrix() {
        System.out.println("-------------------Displaying Theta Matrix--------------------");
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                float c1 = Theta_Matrix[i][j];
                System.out.printf("%-28f", c1);
            }
            System.out.println();
        }
        System.out.println("---------------------------------------------------------------\n\n");
    }

    public static class NewtonRaphson {
        //Jacobian Matrix
        private float[][] Jacobian_Matrix = new float[3][3];

        //matrix declaration to get the inverse of the a matrix
        private float[][] Transpose_Matrix = new float[3][3];
        private float[][] Adjoint_Matrix = new float[3][3];
        private float[][] Inverse_Matrix = new float[3][3];

        private boolean jacobian = false;
        private boolean unknowns = false;
        private int it = 0;
        NewtonRaphson(int iteration) {
            System.out.println("\n<<<<<<<<<<<<<<<<<<<< Power Flow Equation Using Newton Raphson Method >>>>>>>>>>>>>>>>>>>>>>>>\n");
            this.it = iteration;
            for (int i = 0; i < iteration; i++) {
                form_Jacobian_Matrix(i);
                nthIterationMatrix(i);
            }
        }

        NewtonRaphson(int iteration,boolean jacobian) {
            System.out.println("\n<<<<<<<<<<<<<<<<<<<< Power Flow Equation Using Newton Raphson Method >>>>>>>>>>>>>>>>>>>>>>>>\n");
            this.it = iteration;
            this.jacobian = jacobian;
            for (int i = 0; i < iteration; i++) {
                form_Jacobian_Matrix(i);
                nthIterationMatrix(i);
            }
            System.out.println("-----------------------------------------------------------------------\n\n");

        }

        NewtonRaphson(int iteration,boolean jacobian,boolean unknowns) {
            System.out.println("\n<<<<<<<<<<<<<<<<<<<< Power Flow Equation Using Newton Raphson Method >>>>>>>>>>>>>>>>>>>>>>>>\n");
            this.it = iteration;
            this.jacobian = jacobian;
            this.unknowns = unknowns;
            for (int i = 0; i < iteration; i++) {
                form_Jacobian_Matrix(i);
                nthIterationMatrix(i);
            }
            System.out.println("-----------------------------------------------------------------------\n\n");

        }

        private void form_Jacobian_Matrix(int itr) {

            double value = 0.0;

            //J11
            value = Y_Matrix[1][0] * Vol[0] * Vol[1] * HrMath.Sin(Theta_Matrix[1][0] + Angle[0] - Angle[1])
                    + Y_Matrix[1][2] * Vol[2] * Vol[1] * HrMath.Sin(Theta_Matrix[1][2] + Angle[2] - Angle[1]);
            Jacobian_Matrix[0][0] = HrMath.roundOff(value);

            //J12
            value = Y_Matrix[1][2] * Vol[2] * Vol[1] * HrMath.Sin(Theta_Matrix[1][2] + Angle[2] - Angle[1])
                    * -1;
            Jacobian_Matrix[0][1] = HrMath.roundOff(value);

            //J13
            value = Y_Matrix[1][0] * Vol[0] * HrMath.Cos(Theta_Matrix[1][0] + Angle[0] - Angle[1])
                    + (2 * Y_Matrix[1][1] * Vol[1] * HrMath.Cos(Theta_Matrix[1][1]))
                    + Y_Matrix[1][2] * Vol[2] * HrMath.Cos(Theta_Matrix[1][2] + Angle[2] - Angle[1]);
            Jacobian_Matrix[0][2] = HrMath.roundOff(value);

            //J21
            value = Y_Matrix[2][1] * Vol[1] * Vol[2] * HrMath.Sin(Theta_Matrix[2][1] + Angle[1] - Angle[2]) * -1;
            Jacobian_Matrix[1][0] = HrMath.roundOff(value);

            //J22
            value = Y_Matrix[2][0] * Vol[0] * Vol[2] * HrMath.Sin(Theta_Matrix[2][0] + Angle[0] - Angle[2])
                    + Y_Matrix[2][1] * Vol[1] * Vol[2] * HrMath.Sin(Theta_Matrix[2][1] + Angle[1] - Angle[2]);
            Jacobian_Matrix[1][1] = HrMath.roundOff(value);

            //J23
            value = Y_Matrix[2][1] * Vol[2] * HrMath.Cos(Theta_Matrix[2][1] + Angle[1] - Angle[2]);
            Jacobian_Matrix[1][2] = HrMath.roundOff(value);

            //J31
            value = Y_Matrix[1][0] * Vol[0] * Vol[1] * HrMath.Cos(Theta_Matrix[1][0] + Angle[0] - Angle[1])
                    + Y_Matrix[1][2] * Vol[2] * Vol[1] * HrMath.Cos(Theta_Matrix[1][2] + Angle[2] - Angle[1]);
            Jacobian_Matrix[2][1] = HrMath.roundOff(value);

            //J32
            value = (Y_Matrix[1][2] * Vol[2] * Vol[1] * HrMath.Cos(Theta_Matrix[1][2] + Angle[2] - Angle[1])) * -1;
            Jacobian_Matrix[2][1] = HrMath.roundOff(value);

            //J33
            value = (Y_Matrix[1][0] * Vol[0] * HrMath.Sin(Theta_Matrix[1][0] + Angle[0] - Angle[1])
                    + (2 * Y_Matrix[1][1] * Vol[1] * HrMath.Sin(Theta_Matrix[1][1]))
                    + Y_Matrix[1][2] * Vol[2] * HrMath.Sin(Theta_Matrix[1][2] + Angle[2] - Angle[1])) * -1;
            Jacobian_Matrix[2][2] = HrMath.roundOff(value);

            if (jacobian==true) displayJacobianMatrix(itr);

        }

        private float getDeterminantOfMatrix() {
            float determinant, x, y, z;
            x = (Jacobian_Matrix[0][0] * (Jacobian_Matrix[1][1] * Jacobian_Matrix[2][2]
                    - Jacobian_Matrix[1][2] * Jacobian_Matrix[2][1]));
            y = (Jacobian_Matrix[0][1] * (Jacobian_Matrix[1][0] * Jacobian_Matrix[2][2]
                    - Jacobian_Matrix[1][2] * Jacobian_Matrix[2][0]));
            z = (Jacobian_Matrix[0][2] * (Jacobian_Matrix[1][0] * Jacobian_Matrix[2][1]
                    - Jacobian_Matrix[1][1] * Jacobian_Matrix[2][0]));
            determinant = x - y + z;
            return determinant;
        }

        private void getTransposeOfMatrix() {
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    Transpose_Matrix[i][j] = Jacobian_Matrix[j][i];
                }
            }
        }

        private void getAdjointMatrix() {
            //first row
            Adjoint_Matrix[0][0] = Transpose_Matrix[1][1] * Transpose_Matrix[2][2] - Transpose_Matrix[1][2] * Transpose_Matrix[2][1];
            Adjoint_Matrix[0][1] = Transpose_Matrix[1][0] * Transpose_Matrix[2][2] - Transpose_Matrix[1][2] * Transpose_Matrix[2][0];
            Adjoint_Matrix[0][2] = Transpose_Matrix[1][0] * Transpose_Matrix[2][1] - Transpose_Matrix[1][1] * Transpose_Matrix[2][0];
            //second row
            Adjoint_Matrix[1][0] = Transpose_Matrix[0][1] * Transpose_Matrix[2][2] - Transpose_Matrix[0][2] * Transpose_Matrix[2][1];
            Adjoint_Matrix[1][1] = Transpose_Matrix[0][0] * Transpose_Matrix[2][2] - Transpose_Matrix[0][2] * Transpose_Matrix[2][0];
            Adjoint_Matrix[1][2] = Transpose_Matrix[0][0] * Transpose_Matrix[2][1] - Transpose_Matrix[0][1] * Transpose_Matrix[2][0];
            //third row
            Adjoint_Matrix[2][0] = Transpose_Matrix[0][1] * Transpose_Matrix[1][2] - Transpose_Matrix[0][2] * Transpose_Matrix[1][1];
            Adjoint_Matrix[2][1] = Transpose_Matrix[0][0] * Transpose_Matrix[1][2] - Transpose_Matrix[0][2] * Transpose_Matrix[1][0];
            Adjoint_Matrix[2][2] = Transpose_Matrix[0][0] * Transpose_Matrix[1][1] - Transpose_Matrix[0][1] * Transpose_Matrix[1][0];

            //get the cofactors
            Adjoint_Matrix[0][1] = Adjoint_Matrix[0][1] * -1;
            Adjoint_Matrix[1][0] = Adjoint_Matrix[1][0] * -1;
            Adjoint_Matrix[1][2] = Adjoint_Matrix[1][2] * -1;
            Adjoint_Matrix[2][1] = Adjoint_Matrix[2][1] * -1;

        }

        private void getInverseMatrix() {
            float det = getDeterminantOfMatrix();
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    Inverse_Matrix[i][j] = (1 / det) * Adjoint_Matrix[i][j];
                }
            }
        }

        private void nthIterationMatrix(int itr) {

            getTransposeOfMatrix();
            getAdjointMatrix();
            getInverseMatrix();

            float[] resultantProduct = new float[]{0f, 0f, 0f};
            f[0] = fOf("P2");
            f[1] = fOf("P3");
            f[2] = fOf("Q2");
            // Multiply the two matrices
            for (int i = 0; i < 3; i++) {
                for (int k = 0; k < 3; k++)
                    resultantProduct[i] += Inverse_Matrix[i][k] * f[k];
            }

            Angle[1] = Angle[1] - resultantProduct[0];
            Angle[2] = Angle[2] - resultantProduct[1];
            Vol[1] = Vol[1] - resultantProduct[2];

            if (unknowns==true) displayUnknownVal(itr);
        }

        private void displayJacobianMatrix(int itr){

            System.out.println("\n---------------Jacobian matrix for iteration : "+(itr+1)+" ---------------------");

            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    float c1 = Jacobian_Matrix[i][j];
                    System.out.printf("%-28f", c1);
                }
                System.out.println();
            }


        }

        private void displayUnknownVal(int itr){
            System.out.println("\n----------------------------------------------------------------------");
            System.out.println("The values of unknown after iteration " +(itr+1)+" are : ");
            System.out.printf("%-16s%-20f\n","delta 2",Angle[1]);
            System.out.printf("%-16s%-20f\n","delta 3",Angle[2]);
            System.out.printf("%-16s%-20f\n","voltage 2",Vol[1]);
            System.out.println("-----------------------------------------------------------------------");
        }

        public void displayFinalStdValue(){
            System.out.println("------------------------------------"+it+" iteration over-----------------------------------------");
            System.out.println("The final values of unknown in standard units : ");
            System.out.printf("%-16s%-20f\n","delta 2", HrMath.roundOff(Angle[1] * 180 / Math.PI));
            System.out.printf("%-16s%-20f\n","delta 3", HrMath.roundOff(Angle[2] * 180 / Math.PI));
            System.out.printf("%-16s%-20f\n","voltage 2",Vol[1]);
            System.out.println("\n<<<<<<<<<<<<<<<<<<<< Power Flow Equation Using Newton Raphson Method >>>>>>>>>>>>>>>>>>>>>>>>\n");

        }

    }

}

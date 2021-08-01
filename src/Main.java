public class Main {

    public static void main(String[] args) {

        PowerFlow pf = new PowerFlow(10,4,5);
        pf.display_Y_Matrix();
        pf.display_Theta_Matrix();
        PowerFlow.NewtonRaphson newtonRaphson = new PowerFlow.NewtonRaphson(5,true,true);
        newtonRaphson.displayFinalStdValue();

    }

}

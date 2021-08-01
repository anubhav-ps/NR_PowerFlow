import java.math.BigDecimal;
import java.math.RoundingMode;

public class HrMath {

    public static float Sin(float val) {
        return (float) Math.sin(val);
    }

    public static float Cos(float val) {
        return (float) Math.cos(val);
    }

    public static float Tan(float val){
        return (float) Math.tan(val);
    }

    public static float roundOff(Double num) {
        Double toBeTruncated = num;
        Double truncatedDouble = BigDecimal.valueOf(toBeTruncated)
                .setScale(3, RoundingMode.HALF_DOWN)
                .doubleValue();
        return truncatedDouble.floatValue();
    }

}




public class LoggerTest {
    public double bar() {
        return 1.0;
    }

    public double foo() {
        double a = 0;
        for (int i = 0; i < 5; i++) {
            System.out.printf("Print double %f after bar\n", 2.0);
            pt.up.fe.specs.util.SpecsIo.append(new java.io.File("log.txt"), String.format("Logging to a file\n"));
            a += bar();
            pt.up.fe.specs.util.SpecsIo.append(new java.io.File("log.txt"), String.format("Logging again to a file\n"));
            System.out.printf("Printing again\n");
        }
        return a;
    }

    public static void main(String[] args) {
        System.out.printf("Print double %f after foo\n", 2.0);
        pt.up.fe.specs.util.SpecsIo.append(new java.io.File("log.txt"), String.format("Logging to a file\n"));
        new LoggerTest().foo();
        pt.up.fe.specs.util.SpecsIo.append(new java.io.File("log.txt"), String.format("Logging again to a file\n"));
        System.out.printf("Printing again\n");
    }
}


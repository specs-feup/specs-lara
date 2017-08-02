


public class LoggerTest {
    public double bar() {
        return 1.0;
    }

    public double foo() {
        double a = 0;
        for (int i = 0; i < 5; i++) {
            pt.up.fe.specs.util.SpecsIo.append(new java.io.File("log.txt"), String.format("Logging to a file\n"));
            System.out.printf("Print double %f after println\n", 2.0);
            a += bar();
            System.out.printf("Printing again\n");
            pt.up.fe.specs.util.SpecsIo.append(new java.io.File("log.txt"), String.format("Logging again to a file\n"));
        }
        return a;
    }

    public static void main(String[] args) {
        pt.up.fe.specs.util.SpecsIo.append(new java.io.File("log.txt"), String.format("Logging to a file\n"));
        System.out.printf("Print double %f after println\n", 2.0);
        new LoggerTest().foo();
        System.out.printf("Printing again\n");
        pt.up.fe.specs.util.SpecsIo.append(new java.io.File("log.txt"), String.format("Logging again to a file\n"));
    }
}


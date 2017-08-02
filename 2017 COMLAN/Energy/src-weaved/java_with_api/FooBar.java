


public class FooBar {
    public double bar() {
        return 1.0;
    }

    public double foo() {
        double a = 0;
        for (int i = 0; i < 5; i++) {
            double[] kadabra_energy_output_0Before = weaver.kadabra.monitor.energy.EnergyCheck.getEnergyStats();
            a += bar();
            double[] kadabra_energy_output_0After = weaver.kadabra.monitor.energy.EnergyCheck.getEnergyStats();
            double kadabra_energy_output_0 = 0;
            for(int kadabra_energy_output_0Counter = 0; kadabra_energy_output_0Counter < kadabra_energy_output_0Before.length; kadabra_energy_output_0Counter++){
                kadabra_energy_output_0 += kadabra_energy_output_0After[ kadabra_energy_output_0Counter ] - kadabra_energy_output_0Before[ kadabra_energy_output_0Counter ]; // /10?
            }
            //weaver.kadabra.monitor.energy.EnergyCheck.ProfileDealloc();
            System.out.printf("Energy:%f", kadabra_energy_output_0);
        }
        return a;
    }

    public static void main(String[] args) {
        double[] kadabra_energy_output_1Before = weaver.kadabra.monitor.energy.EnergyCheck.getEnergyStats();
        new FooBar().foo();
        double[] kadabra_energy_output_1After = weaver.kadabra.monitor.energy.EnergyCheck.getEnergyStats();
        double kadabra_energy_output_1 = 0;
        for(int kadabra_energy_output_1Counter = 0; kadabra_energy_output_1Counter < kadabra_energy_output_1Before.length; kadabra_energy_output_1Counter++){
            kadabra_energy_output_1 += kadabra_energy_output_1After[ kadabra_energy_output_1Counter ] - kadabra_energy_output_1Before[ kadabra_energy_output_1Counter ]; // /10?
        }
        //weaver.kadabra.monitor.energy.EnergyCheck.ProfileDealloc();
        System.out.printf("Energy:%f", kadabra_energy_output_1);
    }
}


package org.petuum.app.tracking4;

public class Tracking {

	public static void main(String[] args) {
		long start_time = System.currentTimeMillis();
		TrackingEngine reg_engine = new TrackingEngine();
		int clientId = Integer.parseInt(args[0]);
		String configFileName = args[1];

		// For script purpose
//		double lambda = Double.parseDouble(args[2]);
//		double gamma = Double.parseDouble(args[3]);
//		double eta = Double.parseDouble(args[4]);
//		double alpha = Double.parseDouble(args[5]);
//		double beta = Double.parseDouble(args[6]);
//		int K = Integer.parseInt(args[7]);

//		reg_engine.start(clientId, configFileName, lambda, gamma, eta, alpha,
//				beta, K);
		 reg_engine.start(clientId, configFileName);
		System.out.printf("The program took %d milli sec\n",
				System.currentTimeMillis() - start_time);
	}
}

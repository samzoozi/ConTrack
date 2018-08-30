package org.petuum.app.tracking4;

/**
 * 
 * @author yihuaf
 *
 */

import org.apache.commons.configuration.Configuration;
import org.apache.commons.configuration.ConfigurationException;
import org.apache.commons.configuration.PropertiesConfiguration;

public class ConfigurationBuilder {
	public int numCommChannelPerClient;
	public int numClients;
	public int numTables;
	public int numWorkerThreads;
	public String sspMode;
	public int numClockPerIteration;
	public float initStepSize;
	public boolean useStepDec;
	public double stepDec;
	public int numClocksPerEval;
	public int numIterations;
	public int K;
	public double lambda;
	public int staleness;
	public int staleness_B;
	public int cacheSize;
	public String hostfile;
	public boolean useAdaGrad;
	public double stepSize_theta;
	public double stepSize_B;
	public double stepSize_z;
	public double alphagrad;
	public double alphagrad_B;
	public double alphagrad_theta;
	public double alphagrad_z;
	public double alpha;
	public double beta;
	public double eta;
	public double gamma;
	public int dimension;
	public String inputfile;
	public String initB;
	public String initTheta;
	public String B_file;
	public String theta_file;
	public String objective_file;
	public String time_file;
	public String clock_file;
	public String test_file;
	public int num_users;

	public ConfigurationBuilder(String configFileName) {
		LoadConfiguration(configFileName);
	}

	private void LoadConfiguration(String configFileName) {

		try {
			Configuration config = new PropertiesConfiguration(configFileName);

			// Parse configuration from property files.
			numCommChannelPerClient = config.getInt("numCommChannelPerClient",
					2);
			numClients = config.getInt("numClients", 2);
			numTables = config.getInt("numTables", 3);
			numWorkerThreads = config.getInt("numWorkerThreads", 4);
			sspMode = config.getString("sspMode", "SSP");
			numClockPerIteration = config.getInt("numClockPerIteration", 1);
			initStepSize = config.getFloat("initStepSize", 0.05f);
			useStepDec = config.getBoolean("useStepDec", false);
			stepDec = config.getDouble("stepDec", 0.95);
			numClocksPerEval = config.getInt("numClocksPerEval", 10);
			numIterations = config.getInt("numIterations", 1000);
			K = config.getInt("K", 100);
			lambda = config.getDouble("lambda", 1.0);
			staleness = config.getInt("staleness", 0);
			staleness_B = config.getInt("staleness_B", 1);
			cacheSize = config.getInt("cacheSize", 100000);
			hostfile = config.getString("hostfile", "");
			useAdaGrad = config.getBoolean("useAdaGrad", false);
			stepSize_theta = config.getDouble("stepSize_theta", 0.1);
			stepSize_B = config.getDouble("stepSize_B", 0.1);
			stepSize_z = config.getDouble("stepSize_z", 0.1);
			alphagrad = config.getDouble("alphagrad", 0.1);
			alphagrad_B = config.getDouble("alphagrad_B", 0.1);
			alphagrad_theta = config.getDouble("alphagrad_theta", 0.1);
			alphagrad_z = config.getDouble("alphagrad_z", 0.1);
			alpha = config.getDouble("alpha", 1.0);
			beta = config.getDouble("beta", 1.0);
			eta = config.getDouble("eta", 10);
			gamma = config.getDouble("gamma", 10);
			dimension = config.getInt("dimension");
			inputfile = config.getString("inputfile");
			initB = config.getString("initB");
			initTheta = config.getString("initTheta");
			B_file = config.getString("B_file");
			theta_file = config.getString("theta_file");
			objective_file = config.getString("objective_file");
			time_file = config.getString("time_file");
			clock_file = config.getString("clock_file");
			test_file = config.getString("test_file");
			num_users = config.getInt("users");

		} catch (ConfigurationException e) {
			e.printStackTrace();
		}
	}
}

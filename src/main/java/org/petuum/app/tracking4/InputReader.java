package org.petuum.app.tracking4;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

public class InputReader {
	static int number_of_transactions = 0;;

	public static ArrayList<ArrayList<ArrayList<double[]>>> readX(
			String datafile) throws FileNotFoundException, IOException {
		int columns = 0;
		double[] trx;
		ArrayList<double[]> transactions = new ArrayList<double[]>();
		ArrayList<ArrayList<double[]>> users = new ArrayList<ArrayList<double[]>>();
		ArrayList<ArrayList<ArrayList<double[]>>> data = new ArrayList<ArrayList<ArrayList<double[]>>>();
		int globalTrxNumber = 0;
		int globalUserNumber = 0;
		try (BufferedReader br = new BufferedReader(new FileReader(datafile))) {
			String line;
			// if ((line = br.readLine()) != null)
			// if (line.split("=")[0].equals("rows"))
			// rows = Integer.parseInt(line.split("=")[1]);
			// else {
			// System.out
			// .println("rows value for matrix X is not defined");
			// System.exit(1);
			// }
			if ((line = br.readLine()) != null)
				if (line.split("=")[0].equals("columns"))
					columns = Integer.parseInt(line.split("=")[1]);
				else {
					System.out
							.println("columns value for matrix X is not defined");
					System.exit(1);
				}

			trx = new double[columns];
			while ((line = br.readLine()) != null) {
				if (line.trim().equals("time")) {
					transactions.add(trx);
					users.add(transactions);
					data.add(users);
					users = new ArrayList<ArrayList<double[]>>();
					transactions = new ArrayList<double[]>();
					trx = new double[columns];
					globalUserNumber = 0;
					globalTrxNumber = 0;
					continue;
				}
				String[] values = line.split("\\s+");
				int trxNumber = Integer.parseInt(values[1]);
				int userNumber = Integer.parseInt(values[0]);

				if (userNumber == globalUserNumber + 1) {
					transactions.add(trx);
					trx = new double[columns];
					users.add(transactions);
					transactions = new ArrayList<double[]>();
					globalUserNumber = userNumber;
					globalTrxNumber = 0;
					number_of_transactions++;
				}

				if (userNumber == globalUserNumber) {
					if (trxNumber == globalTrxNumber + 1) {
						transactions.add(trx);
						trx = new double[columns];
						globalTrxNumber = trxNumber;
						number_of_transactions++;
					}
					if (trxNumber == globalTrxNumber) {
						trx[Integer.parseInt(values[2])] = Double
								.parseDouble(values[3]);
					} else {
						System.err.println("transactions for user "
								+ globalUserNumber
								+ " are not sequentially ordered!");
						System.exit(1);
					}
				} else {
					System.err.println("user are not sequentially order!");
					System.exit(1);
				}

			}
		}
		return data;
	}

	public static ArrayList<ArrayList<ArrayList<double[]>>> readSensor(
			String datafile, int num_columns) throws FileNotFoundException,
			IOException {
		int columns = num_columns;

		double[] record = new double[columns];

		ArrayList<ArrayList<double[]>> sensors = new ArrayList<ArrayList<double[]>>(
				54);
		ArrayList<ArrayList<ArrayList<double[]>>> data = new ArrayList<ArrayList<ArrayList<double[]>>>();

		for (int i = 0; i < 54; i++) {
			ArrayList<double[]> ar = new ArrayList<double[]>();
			sensors.add(ar);
		}

		try (BufferedReader br = new BufferedReader(new FileReader(datafile))) {
			String line;
			double minute = 0;
			boolean saveData = false;
			while ((line = br.readLine()) != null) {

				String[] values = line.split(",");
				minute = Double.parseDouble(values[0]);
				if ((minute % 60 == 0)) {
					if (saveData) {
						data.add(sensors);
						sensors = new ArrayList<ArrayList<double[]>>(54);
						for (int i = 0; i < 54; i++) {
							ArrayList<double[]> ar = new ArrayList<double[]>();
							sensors.add(ar);
						}
						record = new double[columns];
						saveData = false;
					}
				} else {
					saveData = true;
				}
				int sensor = Integer.parseInt(values[5]);
				if (sensor > 54)
					continue;
				// counter++;
				record = new double[columns];

				for (int i = 0; i < columns; i++) {
					record[i] = Double.parseDouble(values[i + 1]);
				}

				sensors.get(sensor - 1).add(record);

			}

			data.add(sensors);

			// for (int i = 0; i < columns; i++) {
			// avg_values[i] = avg_values[i] / counter;
			// }

			// Iterator<ArrayList<ArrayList<double[]>>> data_itr = data
			// .listIterator();
			// while (data_itr.hasNext()) {
			// ArrayList<ArrayList<double[]>> sensorList = data_itr.next();
			// Iterator<ArrayList<double[]>> sensor_itr = sensorList
			// .listIterator();
			// while (sensor_itr.hasNext()) {
			// ArrayList<double[]> readingList = sensor_itr.next();
			// for (int i = 0; i < readingList.size(); i++) {
			// double[] val = readingList.get(i);
			// for (int j = 0; j < val.length; j++) {
			// val[j] = val[j] / globalMax[j];
			// // (val[j] - globalMin[j])
			// // / (globalMax[j] - globalMin[j]) * 2 - 1;
			// // val[j] = val[j] - avg_values[j];
			// }
			// readingList.set(i, val);
			// }
			// }
			// }

		}

		return data;
	}

	public static ArrayList<ArrayList<ArrayList<double[]>>> readAirData(
			String datafile, int num_columns) throws FileNotFoundException,
			IOException {
		int columns = num_columns;
		final int num_aircrafts = 36;
		double[] record = new double[columns];
		ArrayList<ArrayList<double[]>> aircrafts = new ArrayList<ArrayList<double[]>>(
				num_aircrafts);
		ArrayList<ArrayList<ArrayList<double[]>>> data = new ArrayList<ArrayList<ArrayList<double[]>>>();

		for (int i = 0; i < num_aircrafts; i++) {
			ArrayList<double[]> ar = new ArrayList<double[]>();
			aircrafts.add(ar);
		}

		try (BufferedReader br = new BufferedReader(new FileReader(datafile))) {
			String line;
			int last_epoch = 1;
			while ((line = br.readLine()) != null) {

				String[] values = line.split(",");
				int epoch = Integer.parseInt(values[1]);
				if (epoch == last_epoch + 1) {
					data.add(aircrafts);
					aircrafts = new ArrayList<ArrayList<double[]>>(
							num_aircrafts);
					for (int i = 0; i < num_aircrafts; i++) {
						ArrayList<double[]> ar = new ArrayList<double[]>();
						aircrafts.add(ar);
					}
					record = new double[columns];
				}

				int aircraftId = Integer.parseInt(values[0]);
				if ((aircraftId > num_aircrafts) || (aircraftId < 0))
					continue;
				record = new double[columns];
				for (int i = 0; i < columns; i++) {
					record[i] = Double.parseDouble(values[i + 2]);
				}

				aircrafts.get(aircraftId - 1).add(record);
				last_epoch = epoch;
			}

			data.add(aircrafts);

		}

		return data;
	}

	public static ArrayList<ArrayList<ArrayList<double[]>>> readStockData(
			String datafile, int num_columns) throws FileNotFoundException,
			IOException {
		int columns = num_columns;
		final int num_symbols = 30;
		double[] record = new double[columns];

		ArrayList<ArrayList<ArrayList<double[]>>> data = new ArrayList<ArrayList<ArrayList<double[]>>>();

		try (BufferedReader br = new BufferedReader(new FileReader(datafile))) {
			String line;
			// int last_epoch = 1;

			while ((line = br.readLine()) != null) {

				String[] values = line.split(",");
				int epoch = Integer.parseInt(values[1]);

				if (epoch > data.size()) {
					ArrayList<ArrayList<double[]>> symbols = new ArrayList<ArrayList<double[]>>(
							num_symbols);
					for (int i = 0; i < num_symbols; i++) {
						ArrayList<double[]> ar = new ArrayList<double[]>();
						symbols.add(ar);
					}
					data.add(symbols);
				}
				// if (epoch != last_epoch) {
				// data.add(symbols);
				// aircrafts = new ArrayList<ArrayList<double[]>>(
				// num_aircrafts);
				// for (int i = 0; i < num_aircrafts; i++) {
				// ArrayList<double[]> ar = new ArrayList<double[]>();
				// aircrafts.add(ar);
				// }
				// record = new double[columns];
				// }

				int symbolId = Integer.parseInt(values[0]);
				record = new double[columns];
				for (int i = 0; i < columns; i++) {
					record[i] = Double.parseDouble(values[i + 2]);
				}

				data.get(epoch - 1).get(symbolId - 1).add(record);

				// last_epoch = epoch;
			}

		}

		return data;
	}

	public static ArrayList<ArrayList<ArrayList<double[]>>> readSyntheticData(
			String datafile, int num_columns, int num_users)
			throws FileNotFoundException, IOException {
		int columns = num_columns;
		double[] record = new double[columns];

		ArrayList<ArrayList<double[]>> users = new ArrayList<ArrayList<double[]>>(
				num_users);
		ArrayList<ArrayList<ArrayList<double[]>>> data = new ArrayList<ArrayList<ArrayList<double[]>>>();

		for (int i = 0; i < num_users; i++) {
			ArrayList<double[]> ar = new ArrayList<double[]>();
			users.add(ar);
		}

		try (BufferedReader br = new BufferedReader(new FileReader(datafile))) {
			String line;
			int epoch = 0;
			int last_epoch = 1;
			while ((line = br.readLine()) != null) {

				String[] values = line.split(",");
				epoch = Integer.parseInt(values[0]);
				if (epoch == last_epoch + 1) {
					data.add(users);
					users = new ArrayList<ArrayList<double[]>>(num_users);
					for (int i = 0; i < num_users; i++) {
						ArrayList<double[]> ar = new ArrayList<double[]>();
						users.add(ar);
					}
					record = new double[columns];
				}

				int userId = Integer.parseInt(values[1]);
				if ((userId > num_users) || (userId < 0))
					continue;
				record = new double[columns];
				for (int i = 0; i < columns; i++) {
					record[i] = Double.parseDouble(values[i + 2]);
				}

				users.get(userId - 1).add(record);
				last_epoch = epoch;

			}

			data.add(users);

		}

		return data;

	}

	public static ArrayList<ArrayList<ArrayList<double[]>>> readActivityData(
			String datafile, int num_columns, int num_users)
			throws FileNotFoundException, IOException {
		int columns = num_columns;
		double[] record = new double[columns];

		ArrayList<ArrayList<ArrayList<double[]>>> data = new ArrayList<ArrayList<ArrayList<double[]>>>();

		for (int i = 0; i < 90; i++) {
			ArrayList<ArrayList<double[]>> users = new ArrayList<ArrayList<double[]>>(
					num_users);
			for (int j = 0; j < num_users; j++) {
				ArrayList<double[]> ar = new ArrayList<double[]>();
				users.add(ar);
			}
			data.add(users);
		}

		try (BufferedReader br = new BufferedReader(new FileReader(datafile))) {
			String line;
			int epoch = 0;
			// int last_epoch = 1;
			// int user = 0;
			// int last_user = 1;
			while ((line = br.readLine()) != null) {

				String[] values = line.split(",");
				epoch = Integer.parseInt(values[0]);
				int userId = Integer.parseInt(values[1]);
				// user = Integer.parseInt(values[1]);
				// if (user == last_user + 1) {
				// data.add(users);
				// users = new ArrayList<ArrayList<double[]>>(num_users);
				// for (int i = 0; i < num_users; i++) {
				// ArrayList<double[]> ar = new ArrayList<double[]>();
				// users.add(ar);
				// }
				// record = new double[columns];
				// }

				if ((userId > num_users) || (userId < 0))
					continue;
				record = new double[columns];
				for (int i = 0; i < columns; i++) {
					record[i] = Double.parseDouble(values[i + 2]);
				}
				data.get(epoch - 1).get(userId - 1).add(record);
				// last_epoch = epoch;

			}

		}

		System.out.println(data.get(1).get(0).size());
		return data;

	}

	public static ArrayList<ArrayList<ArrayList<double[]>>> readUsenetData(
			String datafile) throws FileNotFoundException, IOException {
		int columns = 100;
		final int num_users = 1;
		double[] record = new double[columns];
		ArrayList<ArrayList<double[]>> users = new ArrayList<ArrayList<double[]>>(
				num_users);
		ArrayList<ArrayList<ArrayList<double[]>>> data = new ArrayList<ArrayList<ArrayList<double[]>>>();

		for (int i = 0; i < num_users; i++) {
			ArrayList<double[]> ar = new ArrayList<double[]>();
			users.add(ar);
		}

		int userId = 1;
		int line_num = 1;
		try (BufferedReader br = new BufferedReader(new FileReader(datafile))) {
			String line;

			while ((line = br.readLine()) != null) {

				String[] values = line.split(",");

				record = new double[columns];
				for (int i = 0; i < columns; i++) {
					record[i] = Double.parseDouble(values[i]);
				}
				users.get(userId - 1).add(record);

				if (line_num % 300 == 0) {
					data.add(users);
					users = new ArrayList<ArrayList<double[]>>(num_users);
					for (int i = 0; i < num_users; i++) {
						ArrayList<double[]> ar = new ArrayList<double[]>();
						users.add(ar);
					}
				}
				line_num++;
			}

		}

		return data;
	}

	public static void main(String[] args) {
		try {
			// InputReader.readSensor("sensor_1day.txt");
			InputReader.readUsenetData("usenet/usenet1.txt");
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}

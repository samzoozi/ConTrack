package org.petuum.app.tracking4;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.Random;

import org.petuum.ps.row.double_.DoubleRow;
import org.petuum.ps.table.DoubleTable;

public class Initialization {
	private ArrayList<ArrayList<ArrayList<double[]>>> data;
	private int _M, _K, _D, _T, _R;
	private int[][][] recordIndex = new int[_T][_M][];
	final int WINDOW_SIZE = 5;
	private HashSet<Integer> groundtruth_seeds;

	public Initialization(ArrayList<ArrayList<ArrayList<double[]>>> input_data,
			int num_users, int num_concepts, int num_dimensions,
			int num_epochs, int num_records, int[][][] recordIndex) {
		data = input_data;
		_M = num_users;
		_K = num_concepts;
		_D = num_dimensions;
		_T = num_epochs;
		_R = num_records;
		this.recordIndex = recordIndex;
		groundtruth_seeds = new HashSet<Integer>();
	}

	private void generateRandomGroundtruthSeeds(double percentage) {
		Random rand = new Random();

		groundtruth_seeds.add(0);
		int number_of_seeds = (int) (percentage * _T);
		while (groundtruth_seeds.size() < number_of_seeds)
			groundtruth_seeds.add(rand.nextInt(_T));

	}

	public static boolean is_groundtruth_seed(int epoch) {

		final int JUMP_OVER_EPOCHS_STEP = 5;
		if (epoch % JUMP_OVER_EPOCHS_STEP == 0)
			return true;
		else
			return false;

	}

	public double[][][] init_seeds_KmeansPP_singleThread(
			double[][][] grountruth_seeds) {

		double[][][] concepts_init = new double[_T][_K][_D];

		System.out.println("Initializing Seeds");
		Random random = new Random();
		// Initialize c_0 to a random point in [0,1]^K
		ArrayList<double[]> prev_centers = new ArrayList<double[]>(_K);

		ArrayList<ArrayList<double[]>> centersSlidingWindow = new ArrayList<ArrayList<double[]>>(
				WINDOW_SIZE);

		ArrayList<double[]> avg_centers = new ArrayList<double[]>(_K);

		for (int k = 0; k < _K; k++) {
			avg_centers.add(new double[_D]);
		}

		for (int t = 0; t < _T; t++) {
			System.out.println("Epoch:" + t);

			ArrayList<double[]> centers = new ArrayList<double[]>(_K + 1);

			if (is_groundtruth_seed(t)) {
				for (int k = 0; k < _K; k++) {
					double[] concept = grountruth_seeds[t][k];
					centers.add(new double[_D]);
					for (int d = 0; d < _D; d++) {
						centers.get(k)[d] = concept[d];

					}
				}
			} else {

				boolean searchForUser = true;
				int randomUser = -1, randomTrx = -1;
				while (searchForUser) {
					randomUser = random.nextInt(data.get(t).size());
					if (data.get(t).get(randomUser).size() != 0) {
						randomTrx = random.nextInt(data.get(t).get(randomUser)
								.size());
						searchForUser = false;
					}
				}

				centers.add(new double[_D]);
				for (int d = 0; d < _D; d++) {
					// double tmp = random.nextDouble();
					// centers.get(0)[d] = tmp;
					double tmp = data.get(t).get(randomUser).get(randomTrx)[d];
					centers.get(0)[d] = tmp;
				}
				// Initialize c_1, ..., c_K using Kmeans++ strategy
				ArrayList<Double> probabilityVector = new ArrayList<Double>(_R);

				for (int j = 0; j < _R; j++) {
					probabilityVector.add(0.0);
				}

				Double probabilityMass = 0.0;
				for (int k = 1; k < _K; k++) {

					// Prepare list of distances to closest centers
					probabilityMass = KmeansPP_distanceUpdate(centers,
							prev_centers, probabilityVector, t);

					// Pick next center with probability proportional to
					// closestSquareDistances

					double threshold = random.nextDouble() * probabilityMass;
					double runningSum = 0.0;
					int chosenPointT = 0;
					int chosenPointUser = 0;
					int chosenPointTrx = 0;
					datapoints: for (int i = 0; i < _M; i++) {
						int num_transactions = recordIndex[t][i].length;
						for (int j = 0; j < num_transactions; j++) {
							if (threshold <= runningSum
									+ probabilityVector
											.get(recordIndex[t][i][j])) {
								chosenPointT = t;
								chosenPointUser = i;
								chosenPointTrx = j;
								break datapoints;
							}
							runningSum += probabilityVector
									.get(recordIndex[t][i][j]);
						}
					}

					centers.add(new double[_D]);
					for (int d = 0; d < _D; d++) {
						centers.get(k)[d] = data.get(chosenPointT)
								.get(chosenPointUser).get(chosenPointTrx)[d];
						// System.out.print(centers.get(k)[d] + "\t");

					}

				}
			}

			// remove the random initial point from centers
			// centers.remove(0);

			// create average seeds from sliding window

			// align current centers to the centers in previous epoch
			System.out.println("Aligning Seeds");
			for (int k = 0; k < _K; k++) {
				double[] sum = new double[_D];
				for (int i = 0; i < centersSlidingWindow.size(); i++) {
					for (int d = 0; d < _D; d++) {
						sum[d] += centersSlidingWindow.get(i).get(k)[d];
					}
				}
				for (int d = 0; d < _D; d++) {
					avg_centers.get(k)[d] = sum[d]
							/ centersSlidingWindow.size();
				}
			}
			ArrayList<double[]> alignedCenters;
			if (!is_groundtruth_seed(t)) {

				alignedCenters = alignCenters(centers, avg_centers);
				// alignedCenters = alignCenters(centers, prev_centers);
			} else {
				alignedCenters = centers;
			}

			if (centersSlidingWindow.size() < WINDOW_SIZE) {
				centersSlidingWindow.add(alignedCenters);
			} else {
				centersSlidingWindow.remove(0);
				centersSlidingWindow.add(alignedCenters);
			}
			prev_centers = alignedCenters;

			System.out.println("Epoch " + t);

			for (int k = 0; k < _K; k++) {
				// DoubleRowUpdate CRowUpdate = new DenseDoubleRowUpdate(_D);
				for (int d = 0; d < _D; d++) {
					// CRowUpdate.setUpdate(d, alignedCenters.get(k)[d]);
					concepts_init[t][k][d] = alignedCenters.get(k)[d];
					System.out.print(alignedCenters.get(k)[d] + ",");
					// writerSeeds.printf("%f\t", alignedCenters.get(k)[d]);
				}
				// B.batchInc(t * _K + k, CRowUpdate);
				System.out.println();
				// writerSeeds.println();
			}
		}

		return concepts_init;
	}

	public double KmeansPP_distanceUpdate(ArrayList<double[]> centers,
			ArrayList<double[]> prev_centers,
			ArrayList<Double> probabilityVector, int epoch) {
		double probabilityMass = 0.0;
		ArrayList<Double> pointToCentersSquareDistances = new ArrayList<Double>(
				centers.size());
		ArrayList<Double> pointToLastEpochCentersSquareDistance = new ArrayList<Double>(
				prev_centers.size());
		for (int k = 0; k < centers.size(); k++) {
			pointToCentersSquareDistances.add(0.0);

		}
		for (int k = 0; k < prev_centers.size(); k++) {
			pointToLastEpochCentersSquareDistance.add(0.0);
		}

		probabilityMass = 0.0;

		for (int i = 0; i < _M; i++) {
			int num_transactions = recordIndex[epoch][i].length;
			for (int j = 0; j < num_transactions; j++) {
				// Find squared distances from point j to each center
				double minSquareDistanceToCurrentSeeds = computeNearestCenterSqureDistance(
						centers, pointToCentersSquareDistances, i, j, epoch);

				double minSquareDistanceToLastEpochSeeds = 0;

				if (prev_centers.size() != 0)
					minSquareDistanceToLastEpochSeeds = computeNearestCenterSqureDistance(
							prev_centers,
							pointToLastEpochCentersSquareDistance, i, j, epoch);
				double prob = minSquareDistanceToCurrentSeeds
						* Math.exp(-1
								* Math.sqrt(minSquareDistanceToLastEpochSeeds));

				// double prob = minSquareDistanceToCurrentSeeds;
				probabilityVector.set(recordIndex[epoch][i][j], prob);
				probabilityMass += prob;
			}
		}

		return probabilityMass;
	}

	private double computeNearestCenterSqureDistance(
			ArrayList<double[]> centers, ArrayList<Double> squareDistancesList,
			int i, int j, int epoch) {
		for (int k = 0; k < centers.size(); k++) {
			double squareDistance = 0.0;
			for (int d = 0; d < _D; d++) {

				double tmp = data.get(epoch).get(i).get(j)[d]
						- centers.get(k)[d];

				squareDistance += tmp * tmp;
			}

			squareDistancesList.set(k, squareDistance);

		}

		return Collections.min(squareDistancesList);

	}

	private ArrayList<double[]> alignCenters(ArrayList<double[]> centers,
			ArrayList<double[]> avg_centers) {
		ArrayList<double[]> aligned = new ArrayList<double[]>(_K);
		for (int i = 0; i < _K; i++) {
			aligned.add(new double[_D]);
		}
		// ArrayList<Double> pointToCentersSquareDistances = new
		// ArrayList<Double>(
		// centers.size());
		ArrayList<Integer> usedIndexes = new ArrayList<Integer>();

		for (int i = 0; i < centers.size(); i++) {
			int nearestIndex = -1;
			double nearestDistance = Double.POSITIVE_INFINITY;
			for (int k = 0; k < avg_centers.size(); k++) {
				if (usedIndexes.contains(new Integer(k))) {
					continue;
				}
				double[] prevCenter = avg_centers.get(k);
				double squareDistance = 0.0;
				for (int d = 0; d < _D; d++) {
					// System.out.print(data.get(t).get(i).get(j)[d]
					// + "\t");
					double tmp = centers.get(i)[d] - prevCenter[d];
					// double tmp = edges.get(j).data[d]
					// - centers.get(k).get(d);
					squareDistance += tmp * tmp;
				}
				if (squareDistance < nearestDistance) {
					nearestDistance = squareDistance;
					nearestIndex = k;
				}
			}
			usedIndexes.add(new Integer(nearestIndex));
			for (int d = 0; d < _D; d++) {
				aligned.get(nearestIndex)[d] = centers.get(i)[d];
			}
			// aligned.add(centers.get(nearestIndex));
			// centers.remove(nearestIndex);

		}

		return aligned;
	}

	public double[][][][] init_ZLocalTable_closestCommunity(int user_begin,
			int user_end, DoubleTable B, double[][][][] Z) {

		for (int t = 0; t < _T; t++) {
			ArrayList<DoubleRow> CRows = new ArrayList<DoubleRow>(_K);

			for (int k = 0; k < _K; k++) {
				CRows.add(B.get(t * _K + k));
			}

			int userIdx = 0;
			for (int i = user_begin; i < user_end; i++) {
				int num_transactions = data.get(t).get(i).size();
				if (num_transactions != 0) {
					Z[t][userIdx] = new double[num_transactions][_K];
					// ZTemp[t][userIdx] = new double[num_transactions][_K];
					for (int u = 0; u < num_transactions; u++) {
						int closestCom = findClosestCommunity(data.get(t)
								.get(i).get(u), CRows);
						for (int k = 0; k < _K; k++) {
							Z[t][userIdx][u][k] = (k == closestCom ? 0.9
									: 0.1);

						}
						Z[t][userIdx][u] = simplexProjectionReverse(Z[t][userIdx][u]);
					}
				}
				userIdx++;
			}
		}
		return Z;
	}

	public int findClosestCommunity(double[] data, ArrayList<DoubleRow> CRows) {
		double closestSquareDist = Double.POSITIVE_INFINITY;
		int closestCom = 0;
		for (int k = 0; k < _K; k++) {
			double squareDist = 0.0;
			for (int d = 0; d < _D; d++) {
				double tmp = data[d] - CRows.get(k).get(d);
				squareDist += tmp * tmp;
			}
			if (squareDist < closestSquareDist) {
				closestSquareDist = squareDist;
				closestCom = k;
			}
		}
		return closestCom;
	}

	public void init_Theta_followingZ(int user_begin, int user_end,
			double[][][][] Z, double[][][] theta) {

		// Random random = new Random();

		for (int t = 0; t < _T; t++) {
			int userIdx = 0;
			for (int i = user_begin; i < user_end; i++) {
				int num_transactions = data.get(t).get(i).size();
				if (num_transactions != 0) {
					double[] vector = new double[_K];
					for (int u = 0; u < num_transactions; u++) {
						for (int k = 0; k < _K; k++) {
							vector[k] += Z[t][userIdx][u][k];
						}
					}
					for (int k = 0; k < _K; k++) {
						vector[k] /= num_transactions;
						theta[t][userIdx][k] = vector[k];
						// if (global_worker_id == 0)
						// System.out.print(thetaTemp[t][userIdx][k] +
						// "\t");
					}
					// if (global_worker_id == 0)
					// System.out.println();
				}
				userIdx++;
			}
		}

	}

	private double[] simplexProjectionReverse(double[] y) {
		// double[] output = new double[y.length];
		// double[] y = new double[updateVals.length];
		double[] u = new double[_K];
		for (int j = 0; j < y.length; j++) {
			u[j] = y[j];
		}

		Arrays.sort(u);
		int ro = 0;
		int d = u.length - 1;
		for (int j = d; j >= 0; j--) {
			double val = 0;
			for (int i = d; i >= j; i--) {
				val += u[i];
			}
			if (u[j] + ((1 - val) / (d - j + 1)) > 0) {
				ro = d - j;
			}
		}

		double val = 0;
		for (int i = d; i >= d - ro; i--) {
			val += u[i];
		}
		double lambda = (1 - val) / (ro + 1);
		for (int i = 0; i < y.length; i++) {
			if (y[i] + lambda > 0)
				y[i] = y[i] + lambda;
			else
				y[i] = 0;
		}
		return y;
	}

}

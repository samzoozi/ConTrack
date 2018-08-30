package org.petuum.app.tracking4;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Random;
import java.util.Vector;
import java.util.concurrent.BrokenBarrierException;
import java.util.concurrent.CyclicBarrier;

import org.petuum.ps.PsTableGroup;
import org.petuum.ps.common.util.Timer;
import org.petuum.ps.common.util.Utils;
import org.petuum.ps.config.TableConfig;
import org.petuum.ps.config.TableGroupConfig;
import org.petuum.ps.row.double_.DenseDoubleRowUpdate;
import org.petuum.ps.row.double_.DoubleRow;
import org.petuum.ps.row.double_.DoubleRowUpdate;
import org.petuum.ps.table.DoubleTable;

public class TrackingEngine {
	private int clientId;
	private ConfigurationBuilder config;

	private ArrayList<ArrayList<ArrayList<double[]>>> data;
	private int _M, _K, _D, _T, _R;
	int[][][] recordIndex = new int[_T][_M][];

	public class SolveGD implements Runnable {
		private Thread t;
		private CyclicBarrier barrier;
		private int thread_id;
		private DoubleTable B, theta, B_grad, B_sum_grad, objective, CTable;
		private Vector<Float> elapsedTime;
		private Vector<String> clockTimer;
		private double[][][] BTemp, thetaTemp;
		private double[][][][] ZTemp, z_sum;
		private double[][][] z_partial;
		double[][] theta_partial, theta_sum;
		// private double[][] theta_grad, theta_sumgrad;

		boolean printTime = false;

		private int total_num_worker_threads, global_worker_id;
		private PrintWriter writerTheta, writerB, writerObj, writerSeeds;
		float num_x_rows_per_thread;

		public SolveGD(int thread_id, CyclicBarrier process_barrier) {
			this.thread_id = thread_id;
			this.barrier = process_barrier;
			elapsedTime = new Vector<Float>();
			clockTimer = new Vector<String>();
		}

		public void join() {
			try {
				t.join();
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}

		public void start() {
			if (t == null) {
				t = new Thread(this);
				t.start();
			}
		}

		@Override
		public void run() {
			PsTableGroup.registerThread();
			Timer total_timer = new Timer();
			total_num_worker_threads = get_total_num_workers();
			global_worker_id = get_global_worker_id(thread_id);
			if (global_worker_id == 0)
				System.out.println("Start " + total_timer.elapsed() + "s");
			B = PsTableGroup.getDoubleTableOrDie(0);
			B_grad = PsTableGroup.getDoubleTableOrDie(1);
			B_sum_grad = PsTableGroup.getDoubleTableOrDie(2);
			objective = PsTableGroup.getDoubleTableOrDie(3);
			theta = PsTableGroup.getDoubleTableOrDie(4);
			CTable = PsTableGroup.getDoubleTableOrDie(5);

			if (printTime && global_worker_id == 0)
				System.out.println("Generated initial values "
						+ total_timer.elapsed() + "s");

			num_x_rows_per_thread = (float) _M
					/ (float) total_num_worker_threads;
			int data_begin = (int) (global_worker_id * num_x_rows_per_thread);
			int data_end = (int) ((global_worker_id + 1) * num_x_rows_per_thread);

			ZTemp = new double[_T][data_end - data_begin][][];
			thetaTemp = new double[_T][data_end - data_begin][_K];

			try {
				writerSeeds = new PrintWriter("Kmeans++Seeds.txt", "UTF-8");
			} catch (FileNotFoundException | UnsupportedEncodingException e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
			}

			Initialization init = new Initialization(data, _M, _K, _D, _T, _R,
					recordIndex);

			if (global_worker_id == 0) {
				double[][][] gt_seeds = readSeedFile();
				System.out.println("Started Kmeans initialization");
				double[][][] init_seeds = init
						.init_seeds_KmeansPP_singleThread(gt_seeds);
				writeSeedsToTable(init_seeds);

			}
			writerSeeds.close();
			PsTableGroup.globalBarrier();

			init.init_ZLocalTable_closestCommunity(data_begin, data_end, B,
					ZTemp);
			init.init_Theta_followingZ(data_begin, data_end, ZTemp, thetaTemp);

			double[][][] gt_thetas = readThetaGround();
			if (global_worker_id == 1) {
				System.out.println("Here");
			}
			replace_groundtruth_theta(data_begin, data_end, gt_thetas);
			replace_groundtruth_Z(data_begin, data_end, gt_thetas);

			PsTableGroup.globalBarrier();

			try {
				barrier.await();
			} catch (InterruptedException | BrokenBarrierException e) {
				e.printStackTrace();
			}

			try {
				writerTheta = new PrintWriter(config.theta_file, "UTF-8");
				writerB = new PrintWriter(config.B_file, "UTF-8");
				writerObj = new PrintWriter(config.objective_file, "UTF-8");

			} catch (FileNotFoundException | UnsupportedEncodingException e) {
				e.printStackTrace();
			}

			readBThetaTable(data_begin, data_end);

			z_partial = new double[data_end - data_begin][][];
			z_sum = new double[_T][data_end - data_begin][][];
			initZSum(data_begin, data_end);
			theta_partial = new double[data_end - data_begin][_K];
			theta_sum = new double[data_end - data_begin][_K];

			for (int iter = 0; iter < config.numIterations; ++iter) {
				if (thread_id == 0) {
					System.out.println("iteration " + iter + "; time elapsed: "
							+ total_timer.elapsed() + "s");
					if (global_worker_id == 0) {
						elapsedTime.addElement(total_timer.elapsed());

					}
				}

				if (printTime && global_worker_id == 0)
					clockTimer.addElement("Read Z table "
							+ total_timer.elapsed() + "s");

				// printConsoleTheta(iter);
				// printConsoleZ(iter);

				// printConsoleB(iter);

				calculateObjectiveFunc(iter, data_begin, data_end);

				for (int epoch = 0; epoch < _T; epoch++) {

					if (printTime && global_worker_id == 0)
						clockTimer.addElement("Calculated Objective Function "
								+ total_timer.elapsed() + "s");

					computeBGradients(data_begin, data_end, epoch);
					if (printTime && global_worker_id == 0)
						clockTimer
								.addElement("Computed and stored B Gradients "
										+ total_timer.elapsed() + "s");

					theta_partial = new double[data_end - data_begin][_K];
					computeThetaGradients(data_begin, data_end, epoch);
					if (printTime && global_worker_id == 0)
						clockTimer
								.addElement("Computed and stored Theta Gradients "
										+ total_timer.elapsed() + "s");
					// one more clock so everyone gets total values of
					// gradients.

					PsTableGroup.clock();

					if (global_worker_id == 0) {
						updateB(epoch);
						if (printTime)
							clockTimer.addElement("Updated B "
									+ total_timer.elapsed() + "s");
					}

					updateTheta(data_begin, data_end, epoch);
					if (printTime)
						clockTimer.addElement("Updated Theta "
								+ total_timer.elapsed() + "s");

					PsTableGroup.clock();

					readBThetaTable(data_begin, data_end, epoch);

					if (printTime && global_worker_id == 0)
						clockTimer.addElement("Read B,Theta tables "
								+ total_timer.elapsed() + "s");

					computeZGradients(data_begin, data_end, epoch);
					if (printTime && global_worker_id == 0)
						clockTimer.addElement("Computed Z Gradients "
								+ total_timer.elapsed() + "s");

					updateZ(data_begin, data_end, epoch);
					if (printTime && global_worker_id == 0)
						clockTimer.addElement("Updated Z "
								+ total_timer.elapsed() + "s");
					// PsTableGroup.clock();
				}

				if (global_worker_id == 0) {
					writeObjectiveToFile(iter);
					if (printTime)
						clockTimer
								.addElement("Wrote Objective Function to file "
										+ total_timer.elapsed() + "s");
				}

			}

			writeThetaTempToTable(data_begin, data_end);
			PsTableGroup.clock();

			if (global_worker_id == 0) {
				printTheta();
				writeToDisk();
				printB();
			}

			System.out.println("Deregistering threads");
			PsTableGroup.deregisterThread();
			System.out.println("threads deregistered");
			writerTheta.close();
			writerB.close();
			writerObj.close();

			System.out.println("Files closed");

		}

		private void writeThetaTempToTable(int user_begin, int user_end) {
			for (int t = 0; t < _T; t++) {
				int userIdx = 0;
				for (int i = user_begin; i < user_end; i++) {
					DoubleRowUpdate upd = new DenseDoubleRowUpdate(_K);
					for (int j = 0; j < _K; j++) {
						upd.setUpdate(j, thetaTemp[t][userIdx][j]);
					}
					theta.batchInc((t * _M) + i, upd);
					userIdx++;
				}
			}
		}

		private void initZSum(int user_begin, int user_end) {
			for (int t = 0; t < _T; t++) {
				int userIdx = 0;
				for (int m = user_begin; m < user_end; m++) {
					int num_transactions = data.get(t).get(m).size();
					z_sum[t][userIdx] = new double[num_transactions][_K];
					userIdx++;
				}
			}
		}

		private void readBThetaTable(int user_begin, int user_end) {
			BTemp = new double[_T][_K][_D];
			for (int t = 0; t < _T; t++) {
				for (int i = 0; i < _K; i++) {
					DoubleRow row = B.get((t * _K) + i);
					for (int j = 0; j < _D; j++) {
						BTemp[t][i][j] = row.getUnlocked(j);
					}
				}
			}

		}

		private void readBThetaTable(int user_begin, int user_end, int epoch) {

			for (int i = 0; i < _K; i++) {
				DoubleRow row = B.get((epoch * _K) + i);
				for (int j = 0; j < _D; j++) {
					BTemp[epoch][i][j] = row.getUnlocked(j);
				}
			}
		}

		private void calculateObjectiveFunc(int iteration, int user_begin,
				int user_end) {

			double val1_1 = 0;
			for (int t = 0; t < _T; t++) {
				int userIdx = 0;
				for (int i = user_begin; i < user_end; i++) {
					ArrayList<double[]> transactions = data.get(t).get(i);
					for (int j = 0; j < transactions.size(); j++) {
						for (int d = 0; d < _D; d++) {
							double multiply = multiplyMatricesZB(ZTemp,
									userIdx, j, BTemp, d, t);
							val1_1 += (Math.pow(multiply
									- transactions.get(j)[d], 2))
									/ (transactions.size());
						}
					}
					userIdx++;
				}
			}

			val1_1 = val1_1 / (_T * _M * _D);
			if (global_worker_id == 1)
				System.out.println("val1_1=" + val1_1);

			double val1_2 = 0;
			for (int t = 0; t < _T; t++) {
				int userIdx = 0;
				for (int i = user_begin; i < user_end; i++) {
					double[] sum_z = new double[_K];
					int num_transactions = data.get(t).get(i).size();
					if (num_transactions != 0) {
						for (int j = 0; j < num_transactions; j++) {
							for (int k = 0; k < _K; k++) {
								sum_z[k] += ZTemp[t][userIdx][j][k];
							}
						}
						for (int k = 0; k < _K; k++) {
							sum_z[k] /= num_transactions;
							val1_2 += Math.pow(thetaTemp[t][userIdx][k]
									- sum_z[k], 2);
						}
					}
					userIdx++;
				}
			}

			val1_2 = val1_2 / (_T * _M * _K);
			if (global_worker_id == 1)
				System.out.println("val1_2=" + val1_2);

			double val2 = 0;
			for (int t = 0; t < _T - 1; t++) {
				int userIdx = 0;
				for (int i = user_begin; i < user_end; i++) {
					for (int j = 0; j < _K; j++) {
						val2 += Math.pow(thetaTemp[t][userIdx][j]
								- thetaTemp[t + 1][userIdx][j], 2);
					}
					userIdx++;
				}
			}
			val2 = val2 * (config.eta) / (_T * _M * _K);
			if (global_worker_id == 1)
				System.out.println("obj_eta=" + val2);

			double val3 = 0;
			if (global_worker_id == 0) {
				for (int t = 0; t < _T - 1; t++) {
					for (int i = 0; i < _K; i++) {
						for (int j = 0; j < _D; j++) {
							val3 += Math.pow(BTemp[t][i][j]
									- BTemp[t + 1][i][j], 2);
						}
					}
				}
				val3 = val3 * (config.gamma) / (_T * _K * _D);
				if (global_worker_id == 0)
					System.out.println("obj_gamma=" + val3);
			}

			// Dirichlet term for theta
			double val4 = 0;
			for (int t = 0; t < _T; t++) {
				int userIdx = 0;
				for (int i = user_begin; i < user_end; i++) {
					for (int j = 0; j < _K; j++) {
						val4 += Math.log(thetaTemp[t][userIdx][j] + 10e-6);
					}
					userIdx++;
				}
			}
			val4 = val4 * ((1 - config.alpha) / (_T * _M * _K));
			if (global_worker_id == 1)
				System.out.println("obj_alpha=" + val4);
			// Lasso Regularizaion
			double val5 = 0;
			if (global_worker_id == 0) {
				for (int t = 0; t < _T; t++) {
					for (int i = 0; i < _K; i++) {
						for (int j = 0; j < _D; j++) {
							val5 += Math.abs(BTemp[t][i][j]);
						}
					}
				}
				val5 = val5 * (config.lambda / (_T * _K * _D));
				if (global_worker_id == 0)
					System.out.println("obj_lambda=" + val5);
			}

			// Dirichlet term for Z
			double val6 = 0;
			for (int t = 0; t < _T; t++) {
				int userIdx = 0;
				for (int i = user_begin; i < user_end; i++) {
					int num_transactions = data.get(t).get(i).size();
					for (int j = 0; j < num_transactions; j++) {
						for (int k = 0; k < _K; k++) {
							val6 += (Math.log(ZTemp[t][userIdx][j][k] + 10e-6))
									/ num_transactions;
						}
					}
					userIdx++;
				}
			}
			val6 = val6 * ((1 - config.beta) / (_T * _M * _K));
			if (global_worker_id == 1)
				System.out.println("obj_beta=" + val6);

			double objectiveValue = val1_1 + val1_2 + val2 + val3 + val4 + val5
					+ val6;
			DoubleRowUpdate upd = new DenseDoubleRowUpdate(1);
			upd.setUpdate(0, objectiveValue);
			objective.batchInc(iteration, upd);

		}

		private void computeBGradients(int data_begin, int data_end, int epoch) {

			double[][] B_partial = new double[_K][_D];
			int userIdx = 0;
			for (int i = data_begin; i < data_end; i++) {
				ArrayList<double[]> transactions = data.get(epoch).get(i);
				int num_transactions = transactions.size();

				for (int j = 0; j < num_transactions; j++) {
					for (int q = 0; q < _D; q++) {
						double multiply = multiplyMatricesZB(ZTemp, userIdx, j,
								BTemp, q, epoch);
						double value = transactions.get(j)[q];
						double subtraction = multiply - value;
						for (int p = 0; p < _K; p++) {
							double val = (ZTemp[epoch][userIdx][j][p])
									* subtraction;
							B_partial[p][q] += (2 * val)
									/ (_T * _M * _D * num_transactions);
						}
					}

				}
				userIdx++;
			}

			// add time factor
			if ((global_worker_id == 0) && (_T != 1)) {
				for (int p = 0; p < _K; p++) {
					for (int q = 0; q < _D; q++) {
						if (epoch == 0) {
							B_partial[p][q] += 2
									* (config.gamma)
									* (BTemp[epoch][p][q] - BTemp[epoch + 1][p][q])
									/ (_T * _K * _D);
						} else if (epoch == _T - 1) {
							B_partial[p][q] += 2
									* (config.gamma)
									* (BTemp[epoch][p][q] - BTemp[epoch - 1][p][q])
									/ (_T * _K * _D);
						} else {
							B_partial[p][q] += 2
									* (config.gamma)
									* ((2 * BTemp[epoch][p][q])
											- BTemp[epoch - 1][p][q] - BTemp[epoch + 1][p][q])
									/ (_T * _K * _D);
						}

					}
				}

			}

			// aggregate gradients in shared table
			storeBGradients(B_partial, epoch);

		}

		private void computeThetaGradients(int user_begin, int user_end,
				int epoch) {

			int userIdx = 0;
			for (int m = user_begin; m < user_end; m++) {
				int num_transactions = data.get(epoch).get(m).size();
				if (num_transactions != 0) {
					double[] sum_z = new double[_K];
					for (int j = 0; j < num_transactions; j++) {
						for (int k = 0; k < _K; k++) {
							sum_z[k] += ZTemp[epoch][userIdx][j][k];
						}
					}
					for (int k = 0; k < _K; k++) {
						double avg = sum_z[k] / num_transactions;
						theta_partial[userIdx][k] += 2
								* (thetaTemp[epoch][userIdx][k] - avg)
								/ (_T * _M * _K);
					}
				}
				userIdx++;
			}

			// add time factor
			if (_T != 1) {

				userIdx = 0;
				for (int p = user_begin; p < user_end; p++) {
					for (int q = 0; q < _K; q++) {
						if (epoch == 0) {
							theta_partial[userIdx][q] += 2
									* (config.eta)
									* (thetaTemp[epoch][userIdx][q] - thetaTemp[epoch + 1][userIdx][q])
									/ (_T * _M * _K);
						} else if (epoch == _T - 1) {
							theta_partial[userIdx][q] += 2
									* (config.eta)
									* (thetaTemp[epoch][userIdx][q] - thetaTemp[epoch - 1][userIdx][q])
									/ (_T * _M * _K);
						} else {
							theta_partial[userIdx][q] += 2
									* (config.eta)
									* ((2 * thetaTemp[epoch][userIdx][q])
											- thetaTemp[epoch - 1][userIdx][q] - thetaTemp[epoch + 1][userIdx][q])
									/ (_T * _M * _K);
						}

					}
					userIdx++;
				}

			}

			// add Dirichlet prior
			userIdx = 0;
			for (int p = user_begin; p < user_end; p++) {
				for (int q = 0; q < _K; q++) {
					theta_partial[userIdx][q] += (1 - config.alpha)
							/ ((thetaTemp[epoch][userIdx][q] + 10e-6) * _T * _M * _K);
				}
				userIdx++;
			}

			// storeThetaGradients(theta_partial, user_begin, user_end, epoch);

		}

		private void computeZGradients(int user_begin, int user_end, int epoch) {

			// Term 1
			int userIdx = 0;
			for (int m = user_begin; m < user_end; m++) {
				int num_transactions = data.get(epoch).get(m).size();
				if (num_transactions != 0) {
					z_partial[userIdx] = new double[num_transactions][_K];
					// double[] sum_z = new double[_K];
					for (int k = 0; k < _K; k++) {
						double val = 0;
						for (int j = 0; j < num_transactions; j++) {
							val += ZTemp[epoch][userIdx][j][k];
						}
						double gr = 2
								* ((val / (num_transactions * num_transactions)) - (thetaTemp[epoch][userIdx][k] / num_transactions))
								/ (_T * _M * _K);
						for (int j = 0; j < num_transactions; j++) {
							z_partial[userIdx][j][k] += gr;
						}
					}
				}
				userIdx++;
			}

			// Term 2

			userIdx = 0;
			for (int m = user_begin; m < user_end; m++) {
				int num_transactions = data.get(epoch).get(m).size();
				for (int j = 0; j < num_transactions; j++) {
					for (int d = 0; d < _D; d++) {
						double multiply = multiplyMatricesZB(ZTemp, userIdx, j,
								BTemp, d, epoch);
						double value = data.get(epoch).get(m).get(j)[d];
						for (int k = 0; k < _K; k++) {
							z_partial[userIdx][j][k] += ((2 * BTemp[epoch][k][d]) * (multiply - value))
									/ (_T * _M + _D * num_transactions);
						}
					}
				}
				userIdx++;
			}

			// Dirichlet

			userIdx = 0;
			for (int m = user_begin; m < user_end; m++) {
				int num_transactions = data.get(epoch).get(m).size();
				for (int j = 0; j < num_transactions; j++) {
					for (int k = 0; k < _K; k++) {
						z_partial[userIdx][j][k] += (1 - config.beta)
								/ ((ZTemp[epoch][userIdx][j][k] + 10e-6) * _T
										* _M * _K * num_transactions);
					}
				}
				userIdx++;
			}

		}

		private void storeBGradients(double[][] Bpartial, int epoch) {
			// store new gradients

			for (int i = 0; i < _K; i++) {
				DoubleRowUpdate upd = new DenseDoubleRowUpdate(_D);
				for (int j = 0; j < _D; j++) {
					upd.setUpdate(j, Bpartial[i][j]);
				}
				B_grad.batchInc((epoch * _K) + i, upd);
			}

		}

		// private void storeThetaGradients(double[][] theta_partial,
		// int user_begin, int user_end, int epoch) {
		// // store new gradients
		//
		// int userIdx = 0;
		// for (int i = user_begin; i < user_end; i++) {
		// DoubleRowUpdate upd = new DenseDoubleRowUpdate(_K);
		// for (int j = 0; j < _K; j++) {
		// upd.setUpdate(j, theta_partial[userIdx][j]);
		// }
		// theta_grad.batchInc((epoch * _M) + i, upd);
		// userIdx++;
		// }
		//
		// }

		private double multiplyMatricesZB(double[][][][] first,
				int userFromFirst, int recordFromUser, double[][][] second,
				int columnFromSecond, int timePoint) {
			double val = 0;
			try {
				for (int i = 0; i < _K; i++) {
					val += first[timePoint][userFromFirst][recordFromUser][i]
							* second[timePoint][i][columnFromSecond];
				}
			} catch (ArrayIndexOutOfBoundsException e) {
				System.out.println("user=" + userFromFirst + ",record="
						+ recordFromUser);
			}
			return val;
		}

		private void updateB(int epoch) {

			double[][] new_B_gradients = new double[_K][_D];

			for (int i = 0; i < _K; i++) {
				DoubleRow row = B_grad.get((epoch * _K) + i);
				for (int j = 0; j < _D; j++) {
					new_B_gradients[i][j] = row.getUnlocked(j);
				}
			}

			double[][] sum_B_gradients = new double[_K][_D];

			for (int i = 0; i < _K; i++) {
				DoubleRow row = B_sum_grad.get((epoch * _K) + i);
				for (int j = 0; j < _D; j++) {
					sum_B_gradients[i][j] = row.getUnlocked(j);
				}
			}

			double[][] updateVals = new double[_K][_D];

			for (int i = 0; i < _K; i++) {
				for (int j = 0; j < _D; j++) {
					if (config.useAdaGrad) {

						updateVals[i][j] = adaGrad_B(new_B_gradients[i][j],
								sum_B_gradients[i][j]);

					} else {
						updateVals[i][j] = (-1 * config.stepSize_B * new_B_gradients[i][j]);
					}

					// do regularization
					updateVals[i][j] = lassoRegularizarion(BTemp, epoch, i, j,
							updateVals[i][j]);

				}

				// updateVals[t][i] = simplexProjection(updateVals[t][i], t,
				// i, BTemp);
			}

			for (int i = 0; i < _K; i++) {
				DoubleRowUpdate upd = new DenseDoubleRowUpdate(_D);
				for (int j = 0; j < _D; j++) {
					upd.setUpdate(j, updateVals[i][j]);
				}
				B.batchInc((epoch * _K) + i, upd);
			}

			updateBGradientTables(new_B_gradients, epoch);
		}

		private void updateTheta(int user_begin, int user_end, int epoch) {

			double[][] updateVals = new double[user_end - user_begin][_K];

			int userIdx = 0;
			for (int i = user_begin; i < user_end; i++) {
				for (int j = 0; j < _K; j++) {
					if (config.useAdaGrad) {
						// updateVals[userIdx][j] = adaGrad_theta(
						// theta_partial[userIdx][j],
						// theta_sum[userIdx][j]);
						thetaTemp[epoch][userIdx][j] += adaGrad_theta(
								theta_partial[userIdx][j],
								theta_sum[userIdx][j]);

					} else {
						thetaTemp[epoch][userIdx][j] += (-1
								* config.stepSize_theta * theta_partial[userIdx][j]);
						// updateVals[userIdx][j] = (-1 * config.stepSize *
						// theta_partial[userIdx][j]);
					}
					// if (global_worker_id == 0 && userIdx == 0 && epoch == 0)
					// {
					// System.out.printf("%.6f", thetaTemp[epoch][userIdx][j]
					// + updateVals[userIdx][j]);
					// System.out.print("\t");
					// }
				}
				// if (global_worker_id == 0 && userIdx == 0 && epoch == 0)
				// System.out.println();

				// updateVals[userIdx] = simplexProjection(updateVals[userIdx],
				// epoch, userIdx, thetaTemp);
				// thetaTemp[epoch][userIdx] =
				simplexProjection(thetaTemp[epoch][userIdx]);

				userIdx++;
			}


			updateThetaGradientTables(user_begin, user_end, theta_partial,
					epoch);
		}

		private void updateZ(int user_begin, int user_end, int epoch) {

			int userIdx = 0;
			for (int i = user_begin; i < user_end; i++) {
				int num_transactions = data.get(epoch).get(i).size();
				for (int j = 0; j < num_transactions; j++) {
					for (int k = 0; k < _K; k++) {
						if (config.useAdaGrad) {
							double g_new = z_partial[userIdx][j][k];
							z_sum[epoch][userIdx][j][k] += Math.pow(g_new, 2);
							ZTemp[epoch][userIdx][j][k] += -1
									* config.alphagrad_z
									* g_new
									/ (Math.sqrt(z_sum[epoch][userIdx][j][k] + 10e-6));
						} else {
							ZTemp[epoch][userIdx][j][k] += (-1
									* config.stepSize_z * z_partial[userIdx][j][k]);
						}
					}
					ZTemp[epoch][userIdx][j] = simplexProjection(ZTemp[epoch][userIdx][j]);
				}
				userIdx++;
			}

		}

		private double lassoRegularizarion(double[][][] oldValues, int t,
				int i, int j, double updateValue) {
			double newValue = oldValues[t][i][j] + updateValue;
			double coef = config.lambda / (_T * _K * _D);
			if (newValue > coef)
				return updateValue - coef;
			else if (newValue < -1 * coef)
				// return -1 * (oldValues[t][i][j]);
				return updateValue + coef;
			else
				return -1 * (oldValues[t][i][j]);
		}

		private double nonNegativeProjection(double[][][] oldValues, int t,
				int i, int j, double updateValue) {
			double newValue = oldValues[t][i][j] + updateValue;
			if (newValue > config.lambda)
				return updateValue - config.lambda;
			else
				return -1 * (oldValues[t][i][j]);
		}

		private double[] simplexProjection(double[] updateVals, int _t, int _i,
				double[][][] oldValue) {
			double[] output = new double[updateVals.length];
			double[] y = new double[updateVals.length];
			Double[] u = new Double[updateVals.length];
			for (int j = 0; j < updateVals.length; j++) {
				y[j] = oldValue[_t][_i][j] + updateVals[j];
				u[j] = new Double(y[j]);
			}
			Arrays.sort(u, Collections.reverseOrder());
			int ro = 0;
			for (int j = 0; j < u.length; j++) {
				double val = 0;
				for (int i = 0; i <= j; i++) {
					val += u[i];
				}
				if (u[j] + ((1 - val) / (j + 1)) > 0)
					ro = j;
			}

			double val = 0;
			for (int i = 0; i <= ro; i++) {
				val += u[i];
			}
			double lambda = (1 - val) / (ro + 1);
			for (int i = 0; i < updateVals.length; i++) {
				if (y[i] + lambda > 0)
					output[i] = updateVals[i] + lambda;
				else
					output[i] = -1 * oldValue[_t][_i][i];
			}
			return output;
		}

		private double[] simplexProjection(double[] updateVals, int _t, int _i,
				int _j, double[][][][] oldValue) {

			double[] output = new double[updateVals.length];
			double[] y = new double[updateVals.length];
			Double[] u = new Double[updateVals.length];
			for (int k = 0; k < updateVals.length; k++) {
				y[k] = oldValue[_t][_i][_j][k] + updateVals[k];
				u[k] = new Double(y[k]);
			}
			Arrays.sort(u, Collections.reverseOrder());
			int ro = 0;
			for (int j = 0; j < u.length; j++) {
				double val = 0;
				for (int i = 0; i <= j; i++) {
					val += u[i];
				}
				if (u[j] + ((1 - val) / (j + 1)) > 0)
					ro = j;
			}

			double val = 0;
			for (int i = 0; i <= ro; i++) {
				val += u[i];
			}
			double lambda = (1 - val) / (ro + 1);
			for (int i = 0; i < updateVals.length; i++) {
				if (y[i] + lambda > 0)
					output[i] = updateVals[i] + lambda;
				else
					output[i] = -1 * oldValue[_t][_i][_j][i];
			}
			return output;
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

		private double[] simplexProjection(double[] y) {
			// double[] output = new double[y.length];
			// double[] y = new double[updateVals.length];
			Double[] u = new Double[_K];
			for (int j = 0; j < y.length; j++) {
				u[j] = y[j];
			}

			Arrays.sort(u, Collections.reverseOrder());
			int ro = 0;
			for (int j = 0; j < u.length; j++) {
				double val = 0;
				for (int i = 0; i <= j; i++) {
					val += u[i];
				}
				if (u[j] + ((1 - val) / (j + 1)) > 0)
					ro = j;
			}

			double val = 0;
			for (int i = 0; i <= ro; i++) {
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

		private double adaGrad_B(double g_new, double g_sum) {

			double val = -1 * config.alphagrad_B * g_new
					/ (Math.sqrt(g_sum + Math.pow(g_new, 2) + 10e-6));

			return val;

		}

		private double adaGrad_theta(double g_new, double g_sum) {

			double val = -1 * config.alphagrad_theta * g_new
					/ (Math.sqrt(g_sum + Math.pow(g_new, 2) + 10e-6));
			return val;

		}

		private void updateBGradientTables(double[][] new_B_gradients, int epoch) {
			// move new gradients to old gradients

			for (int i = 0; i < _K; i++) {
				DoubleRowUpdate upd = new DenseDoubleRowUpdate(_D);
				for (int j = 0; j < _D; j++) {
					upd.setUpdate(j, Math.pow(new_B_gradients[i][j], 2));
				}
				B_sum_grad.batchInc((epoch * _K) + i, upd);
			}

			// set new gradients to zero

			for (int i = 0; i < _K; i++) {
				DoubleRowUpdate upd = new DenseDoubleRowUpdate(_D);
				for (int j = 0; j < _D; j++) {
					upd.setUpdate(j, -1 * new_B_gradients[i][j]);
				}
				B_grad.batchInc((epoch * _K) + i, upd);
			}

		}

		private void updateThetaGradientTables(int user_begin, int user_end,
				double[][] new_theta_gradients, int epoch) {
			// move new gradients to old gradients

			int userIdx = 0;
			for (int i = user_begin; i < user_end; i++) {
				for (int j = 0; j < _K; j++) {
					theta_sum[userIdx][j] += Math.pow(
							new_theta_gradients[userIdx][j], 2);
				}
				userIdx++;
			}

			// set new gradients to zero
			theta_partial = new double[user_end - user_begin][_K];

		}

		private void writeObjectiveToFile(int iter) {
			DoubleRow row = objective.get(iter);
			double val = row.getUnlocked(0);
			writerObj.println(iter + "\t" + val);
			System.out
					.println("Iteration " + iter + " objective value: " + val);

		}

		private void writeToDisk() {
			PrintWriter writerStat = null;
			System.out.println("Writing stats");
			try {
				writerStat = new PrintWriter(config.time_file, "UTF-8");
				for (int i = 0; i < config.numIterations; i++) {
					writerStat.print(i + "\t");
					// DoubleRow row = objective.get(i);
					// writerStat
					// .print(Double.toString(row.getUnlocked(0)) + "\t");
					writerStat.println(Double.toString(elapsedTime.get(i)));
				}

			} catch (FileNotFoundException e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
			} catch (UnsupportedEncodingException e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
			}
			writerStat.close();

		}

		private void printTheta() {

			System.out.println("Writing Theta");

			for (int t = 0; t < _T; t++) {
				// writerTheta.println("t=" + t);

				for (int i = 0; i < _M; i++) {
					DoubleRow r = theta.get((_M * t) + i);
					for (int j = 0; j < _K; j++) {
						writerTheta.printf("%.6f\t", r.getUnlocked(j));
					}
					writerTheta.println();

				}
			}

		}

		private void printB() {

			System.out.println("Writing B");
			for (int t = 0; t < _T; t++) {
				// writerB.println("t=" + t);
				for (int i = 0; i < _K; i++) {
					DoubleRow r = B.get((_K * t) + i);
					for (int j = 0; j < _D; j++) {
						writerB.printf("%f\t", r.getUnlocked(j));
					}
					writerB.println();
				}
			}

		}

		private void printConsoleTheta(int iter) {
			if (global_worker_id == 1) {
				System.out.println("Theta");
				// for (int k = 0; k < _K; k++) {
				// System.out.printf("%.6f\t", (float) (k + 1));
				// }
				// System.out.println();
				// for (int i = 0; i < _T; i++) {
				for (int k = 0; k < _K; k++) {
					System.out.printf("%.6f\t", thetaTemp[0][2][k]);
				}
				// System.out.println();
				// }
				System.out.println();
			}

		}

		private void printConsoleZ(int iter) {
			if (global_worker_id == 1) {
				System.out.println("Z");

				// int num_trans = data.get(0).get(2).size();
				for (int u = 0; u < 10; u++) {
					for (int k = 0; k < _K; k++) {
						System.out.printf("%.6f\t", ZTemp[0][2][u][k]);
					}
					System.out.println();
				}
				System.out.println();

			}
		}

		private void printConsoleB(int iter) {
			if (global_worker_id == 0) {
				System.out.println("B");
				for (int t = 0; t < 5; t++) {
					for (int d = 0; d < _D; d++) {
						System.out.printf("%.6f\t", BTemp[t][0][d]);
					}
					System.out.println();
				}

				// for (int d = 0; d < _D; d++) {
				// System.out.printf("%.6f\t", BTemp[1][0][d]);
				// }
				// System.out.println();
			}
		}

		private int get_global_worker_id(int thread_id2) {
			return (clientId * config.numWorkerThreads) + thread_id;
		}

		private int get_total_num_workers() {
			return config.numClients * config.numWorkerThreads;
		}

		protected void generateTheta(int user_begin, int user_end) {
			double[][][] thetaGround = readThetaGround();
			double[] vector = new double[_K];
			double newVal = 0;
			for (int t = 0; t < _T / 5; t++) {
				int userIdx = 0;
				for (int i = user_begin; i < user_end; i++) {
					for (int j = 0; j < _K; j++) {
						// vector[j] = (new Random().nextGaussian() * 0.05 +
						// thetaGround[t][i][j]);

						// vector[j] = (new Random().nextDouble());
						// while (vector[j] == 0)
						// vector[j] = (new Random().nextDouble());

						newVal = thetaGround[t][i][j];
						thetaTemp[t][userIdx][j] = newVal;

					}

					// double[] projected = simplexProjectionReverse(vector);
					// for (int k = 0; k < _K; k++) {
					// thetaTemp[t][userIdx][k] = projected[k];
					// }
					userIdx++;
				}

			}

		}

		protected void generateZ(int user_begin, int user_end) {

			double[][][] thetaGround = readThetaGround();
			double[] vector = new double[_K];

			for (int t = 0; t < _T / 5; t++) {
				int userIdx = 0;
				for (int i = user_begin; i < user_end; i++) {
					int num_transactions = data.get(t).get(i).size();
					ZTemp[t][userIdx] = new double[num_transactions][_K];
					for (int u = 0; u < num_transactions; u++) {
						for (int j = 0; j < _K; j++) {
							vector[j] = (thetaGround[t][i][j]);
						}
						double[] projected = simplexProjectionReverse(vector);
						for (int k = 0; k < _K; k++) {
							ZTemp[t][userIdx][u][k] = projected[k];
						}

					}
					userIdx++;
				}
			}

		}

		private double[][][] readThetaGround() {
			double[][][] theta = new double[_T][_M][_K];
			try (BufferedReader br = new BufferedReader(new FileReader(
					config.initTheta))) {
				String line;
				try {
					for (int t = 0; t < _T; t++) {
						for (int i = 0; i < _M; i++) {
							line = br.readLine();
							String[] points = line.split("\t");
							for (int j = 0; j < _K; j++) {
								theta[t][i][j] = Double.parseDouble(points[j]);
							}
						}
					}
				} catch (IOException e) {
					e.printStackTrace();
				}
			} catch (FileNotFoundException e1) {
				e1.printStackTrace();
				System.exit(1);
			} catch (IOException e1) {
				e1.printStackTrace();
			}
			return theta;

		}

		private double[][][] readSeedFile() {
			double[][][] seeds = new double[_T][_K][_D];
			try (BufferedReader br = new BufferedReader(new FileReader(
					config.initB))) {
				String line;
				for (int t = 0; t < _T; t++) {
					for (int k = 0; k < _K; k++) {
						line = br.readLine();
						String[] points = line.split(",");
						for (int d = 0; d < _D; d++) {
							seeds[t][k][d] = Double.parseDouble(points[d]);
						}
					}
				}

			} catch (FileNotFoundException e1) {
				e1.printStackTrace();
				System.exit(1);
			} catch (IOException e1) {
				e1.printStackTrace();
			}
			return seeds;
		}

		private void replace_groundtruth_theta(int user_begin, int user_end,
				double[][][] thetaGround) {

			for (int t = 0; t < _T; t++) {
				if (Initialization.is_groundtruth_seed(t)) {
					int userIdx = 0;
					for (int i = user_begin; i < user_end; i++) {
						for (int j = 0; j < _K; j++) {
							thetaTemp[t][userIdx][j] = thetaGround[t][i][j];
						}

						userIdx++;
					}
				}
			}
		}

		private void replace_groundtruth_Z(int user_begin, int user_end,
				double[][][] thetaGround) {

			for (int t = 0; t < _T; t++) {
				if (Initialization.is_groundtruth_seed(t)) {
					int userIdx = 0;
					for (int i = user_begin; i < user_end; i++) {
						int num_transactions = data.get(t).get(i).size();
						// ZTemp[t][userIdx] = new double[num_transactions][_K];
						for (int u = 0; u < num_transactions; u++) {
							for (int j = 0; j < _K; j++) {
								ZTemp[t][userIdx][u][j] = (thetaGround[t][i][j]);
							}

						}
						userIdx++;
					}
				}
			}
		}

		private void writeSeedsToTable(double[][][] seeds) {
			for (int t = 0; t < _T; t++) {
				for (int k = 0; k < _K; k++) {
					DoubleRowUpdate CRowUpdate = new DenseDoubleRowUpdate(_D);
					for (int d = 0; d < _D; d++) {
						CRowUpdate.setUpdate(d, seeds[t][k][d]);
					}
					B.batchInc(t * _K + k, CRowUpdate);
					// System.out.println();
					// writerSeeds.println();
				}
			}
		}

	}

	// public void start(int clientId, String configFileName, double lambda,
	// double gamma, double eta, double alpha, double beta, int K) {
	public void start(int clientId, String configFileName) {

		System.out.println("Starting client id=" + clientId);
		this.clientId = clientId;
		System.out.println("Loading configuration from " + configFileName);
		this.config = new ConfigurationBuilder(configFileName);
		TableGroupConfig table_group_config = new TableGroupConfig();
		table_group_config
				.setNumCommChannelsPerClient(config.numCommChannelPerClient)
				.setNumTotalClients(config.numClients).setNumTables(6)
				.setNumLocalAppThreads(config.numWorkerThreads + 1)
				.setClientId(clientId);
		Utils.getHostInfos(config.hostfile, table_group_config.getHostMap());

		System.out.println("PsTableGroup init!");
		PsTableGroup.init(table_group_config, false);

		System.out.println("PsTableGroup Initialized");

		System.out.println("Reading data");
		Timer load_timer = new Timer();

		try {
			_D = config.dimension;
			// data = InputReader.readSensor(config.inputfile, _D);
			// data = InputReader.readAirData(
			// "air_data/air_normalized_epoch11-44_cont_features.txt", _D);
			// data = InputReader.readStockData(config.inputfile, _D);
			// data = InputReader.readSyntheticData(config.inputfile, _D,
			// config.num_users);
			data = InputReader.readActivityData(config.inputfile, _D,
					config.num_users);
			_T = data.size();
			_M = config.num_users;
			_K = config.K;
		} catch (IOException e) {
			e.printStackTrace();
		}

		_R = 0;
		int idx = 0;
		recordIndex = new int[_T][_M][];
		for (int t = 0; t < data.size(); t++) {
			for (int i = 0; i < _M; i++) {
				int num_transactions = data.get(t).get(i).size();
				_R += num_transactions;
				recordIndex[t][i] = new int[num_transactions];
				for (int j = 0; j < num_transactions; j++) {
					recordIndex[t][i][j] = idx;
					idx++;
				}
			}
		}

		System.out.println("Epochs: " + _T);
		System.out.println("Users: " + _M);
		System.out.println("Records: " + _R);
		// for script purspose
		// config.lambda = lambda;
		// config.gamma = gamma;
		// config.eta = eta;
		// config.alpha = alpha;
		// config.beta = beta;
		// _K = K;

		System.out.println("lambda =" + config.lambda);
		System.out.println("gamma =" + config.gamma);
		System.out.println("eta =" + config.eta);
		System.out.println("alpha =" + config.alpha);
		System.out.println("K =" + _K);

		System.out.println("Loaded data in " + load_timer.elapsed() + "s");

		System.out.println("Reading data done");

		TableConfig table_config_staleness_0 = new TableConfig();

		table_config_staleness_0.setStaleness(config.staleness)
				.setProcessCacheCapacity(config.cacheSize)
				.setNoOplogReplay(false);

		TableConfig table_config_staleness_B = new TableConfig();

		table_config_staleness_B.setStaleness(config.staleness_B)
				.setProcessCacheCapacity(config.cacheSize)
				.setNoOplogReplay(false);

		System.out.println("Creating tables");

		PsTableGroup.createDenseDoubleTable(0, _D, table_config_staleness_B); // B

		PsTableGroup.createDenseDoubleTable(1, _D, table_config_staleness_B); // reg_B

		PsTableGroup.createDenseDoubleTable(2, _D, table_config_staleness_B); // sum_reg_B

		PsTableGroup.createDenseDoubleTable(3, 1, table_config_staleness_0); // objective

		PsTableGroup.createDenseDoubleTable(4, _K, table_config_staleness_0); // Theta

		PsTableGroup.createDenseDoubleTable(5, _D, table_config_staleness_0); // Initialization
																				// points
		System.out.println("Creating tables done");

		// table_config.setProcessCacheCapacity(10000);
		// PsTableGroup.createDenseDoubleTable(2, 6, table_config);

		PsTableGroup.createTableDone();

		SolveGD[] threads = new SolveGD[config.numWorkerThreads];
		CyclicBarrier barrier = new CyclicBarrier(config.numWorkerThreads);
		System.out.println("Creating worker threads!");
		// run
		for (int i = 0; i < config.numWorkerThreads; i++) {
			threads[i] = new SolveGD(i, barrier);
			threads[i].start();
		}
		// join
		for (int i = 0; i < config.numWorkerThreads; i++) {
			threads[i].join();
		}

		System.out.println("System shutting down!");

		PsTableGroup.shutdown();

		return;
	}
}

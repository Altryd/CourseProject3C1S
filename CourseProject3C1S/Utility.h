#pragma warning(disable : 4996)
#pragma once
#include <array>
#include <string>
#include <vector>
#include <openssl/sha.h>
std::array<double, 4> sha256_ret_keys(const std::string& str)
{
	unsigned char hash[SHA256_DIGEST_LENGTH];
	SHA256_CTX sha256;
	SHA256_Init(&sha256);
	SHA256_Update(&sha256, str.c_str(), str.size());
	SHA256_Final(hash, &sha256);
	const short* reinterpreted_const_short = reinterpret_cast<const short*>(hash);
	const int* reinterpreted_const_int = reinterpret_cast<const int*>(hash);
	std::array<double, 4> keys;
	double int_part;
	for (size_t i = 0; i < 6; i += 2)
	{
		double to_array = modf((double)reinterpreted_const_short[i] / reinterpreted_const_short[i + 1], &int_part);
		keys[i / 2] = to_array;
	}
	keys[3] = reinterpreted_const_int[3] % 1000 + 1000;
	return keys;
}


std::array<std::vector<double>, 3> get_chaotic_sequences(const std::array<double, 4>& initial_keys, const int N)
{
	std::array<std::vector<double>, 3> C_sequences;
	std::array<double, 4> values_past;
	values_past = initial_keys;
	std::copy(initial_keys.begin(), initial_keys.begin() + 3, values_past.begin());
	double h = 0.001;  //step for Chen system
	double a = 35, b = 3, c = 28;  // parameters for chen system
	std::array<double, 4> values_present = values_past;
	for (size_t i = 0; i < static_cast<size_t>(initial_keys.back()); ++i)
	{
		// xt = x0 + h*(35*(y0-x0))
		// yt = y0+h*(-7*x0-x0*z0+28*y0)
		//zt = z0 + h * (x0 * y0 - 3 * z0)
		values_present[0] = values_past[0] + h * (a * (values_past[1] - values_past[0]));
		values_present[1] = values_past[1] + h * ((c - a) * values_past[0] - values_past[0] * values_past[2] + 28 * values_past[1]);
		values_present[2] = values_past[2] + h * (values_past[0] * values_past[1] - 3 * values_past[2]);
		std::swap(values_past, values_present);
	}
	size_t iteration_count = static_cast<size_t>(N);
	iteration_count *= iteration_count;
	for (size_t i = 0; i < iteration_count * 3; ++i)
	{
		values_present[0] = values_past[0] + h * (a * (values_past[1] - values_past[0]));
		values_present[1] = values_past[1] + h * ((c - a) * values_past[0] - values_past[0] * values_past[2] + 28 * values_past[1]);
		values_present[2] = values_past[2] + h * (values_past[0] * values_past[1] - 3 * values_past[2]);
		std::swap(values_past, values_present);
		C_sequences[0].push_back(values_present[0]);   //C1
		C_sequences[1].push_back(values_present[1]);   //C2
		C_sequences[2].push_back(values_present[2]);   //C3
	}

	return C_sequences;
}

template <typename T>
std::vector<size_t> sort_ranking(const std::vector<T>& vec)
{
	std::vector<T> copy(vec);
	std::vector<size_t> result(copy.size());
	std::sort(copy.begin(), copy.end());
	for (size_t i = 0; i < vec.size(); ++i)
	{
		auto it = std::lower_bound(copy.begin(), copy.end(), vec[i]);
		auto distance = std::distance(copy.begin(), it);
		result[i] = distance;
	}
	return result;
}


void fill_chaotic_sorting_martrix(int* data, const int start_col, const int start_row,
	const int initial_cols, const int size, const int cols, const int start_val, const int fsm[4])
{
	if (size == 1)
	{
		data[start_row * initial_cols + start_col] = start_val;
		return;
	}
	fill_chaotic_sorting_martrix(data, start_col, start_row, initial_cols, size / 4, cols / 2, start_val + (fsm[0] - 1) * size / 4, fsm);
	fill_chaotic_sorting_martrix(data, start_col + cols / 2, start_row, initial_cols, size / 4, cols / 2, start_val + (fsm[1] - 1) * size / 4, fsm);
	fill_chaotic_sorting_martrix(data, start_col, start_row + cols / 2, initial_cols, size / 4, cols / 2, start_val + (fsm[2] - 1) * size / 4, fsm);
	fill_chaotic_sorting_martrix(data, start_col + cols / 2, start_row + cols / 2, initial_cols, size / 4, cols / 2, start_val + (fsm[3] - 1) * size / 4, fsm);
}
#pragma once
#include <vector>
#include <openssl/sha.h>
#include <opencv2/opencv.hpp>
#include <filesystem>
#include "Utility.h"
class ImageDecryptor
{
private:
	cv::Mat image;
	std::string keys_path;
public:
	ImageDecryptor(const cv::Mat& image, const std::string& keys_path) : image(image), keys_path(keys_path) {}
	bool decrypt_image() const
	{
		if (image.rows % 2 != 0 || image.cols != image.rows)
		{
			throw std::runtime_error("Image Not Found!");
		}
		int dimensions_needed = image.rows;
		double* keys = new double[4];
		FILE* fptr_read = fopen(keys_path.c_str(), "rb");
		fread(keys, sizeof(double), 4, fptr_read);
		fclose(fptr_read);
		std::array<double, 4> keys_array({ keys[0], keys[1], keys[2], keys[3] });
		delete[] keys;
		auto sequences = get_chaotic_sequences(keys_array, dimensions_needed);
		std::vector<std::vector<size_t>> three_matrixes(3);
		auto C1 = sequences[0];
		sequences[0].clear();
		for (size_t i = 0; i < 3; ++i)
		{
			size_t begin = static_cast<size_t>(dimensions_needed) * static_cast<size_t>(dimensions_needed) * i;
			std::vector<double> test(C1.begin() + begin, C1.begin() + begin + static_cast<size_t>(dimensions_needed) * dimensions_needed);
			three_matrixes[i] = sort_ranking(test);
		}


		int fsm[4] = { 4, 2, 1, 3 };
		int* fsm_An = new int[static_cast<size_t>(dimensions_needed) * static_cast<size_t>(dimensions_needed)];
		fill_chaotic_sorting_martrix(fsm_An, 0, 0, dimensions_needed, dimensions_needed * dimensions_needed, dimensions_needed, 1, fsm);

		cv::Mat channels[3];
		cv::Mat decrypted_channels[3];
		cv::Mat will_be_decrypted = image.clone();

		cv::split(image, channels);
		cv::split(will_be_decrypted, decrypted_channels);
		for (size_t channel = 0; channel < 3; ++channel)
		{

			for (size_t i = 0; i < dimensions_needed; ++i)
			{
				for (size_t j = 0; j < dimensions_needed; ++j)
				{
					auto row_s1 = three_matrixes[channel][i * dimensions_needed + j] / dimensions_needed;
					auto col_s1 = three_matrixes[channel][i * dimensions_needed + j] % dimensions_needed;
					auto row_fsm = (fsm_An[i * dimensions_needed + j] - 1) / dimensions_needed;
					auto col_fsm = (fsm_An[i * dimensions_needed + j] - 1) % dimensions_needed;
					decrypted_channels[channel].at<uchar>(row_fsm, col_fsm) = channels[channel].at<uchar>(row_s1, col_s1);
				}
			}
		}
		cv::Mat test = cv::Mat::zeros(image.rows, image.cols, CV_8UC3);
		cv::merge(decrypted_channels, 3, test);
		cv::imwrite("decrypted.png", test);
		return true;
	}
};
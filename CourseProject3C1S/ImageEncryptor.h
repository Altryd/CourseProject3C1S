#include <opencv2/opencv.hpp>
#include <string>
#include "Utility.h"
#pragma once

class ImageEncryptor
{
private:
	cv::Mat image;
	std::string string_to_hash;
public:
	ImageEncryptor(const cv::Mat& image, const std::string& str) : image(image), string_to_hash(str) {}
	bool encrypt_image() const
	{
		if (image.empty())
		{
			throw std::runtime_error("Image Not Found!");
		}
		int N = 0;
		int dimensions_needed = 1;
		if (image.rows > image.cols) N = image.rows;
		else N = image.cols;
		while (dimensions_needed < N)
		{
			dimensions_needed *= 2;
		}
		cv::Mat image_to_resize = cv::Mat::zeros(dimensions_needed, dimensions_needed, image.type());
		for (int i = 0; i < image.rows; ++i)
		{
			for (int j = 0; j < image.cols; ++j)
			{
				image_to_resize.at<cv::Vec3b>(i, j) = cv::Vec3b(image.at<cv::Vec3b>(i, j));
			}
		}
		std::array<double, 4> keys_array = sha256_ret_keys(string_to_hash);
		auto sequences = get_chaotic_sequences(keys_array, dimensions_needed);
		std::vector<std::vector<size_t>> three_matrixes(3);
		auto C1 = sequences[0];
		sequences[0].clear();
		for (size_t i = 0; i < 3; ++i)
		{
			size_t begin = static_cast<size_t>(dimensions_needed) * static_cast<size_t>(dimensions_needed) * i;
			std::vector<double> test(C1.begin() + begin,
				C1.begin() + begin + static_cast<size_t>(dimensions_needed) * static_cast<size_t>(dimensions_needed));
			three_matrixes[i] = sort_ranking(test);
		}

		int fsm[4] = { 4, 2, 1, 3 };
		int* fsm_An = new int[static_cast<size_t>(dimensions_needed) * static_cast<size_t>(dimensions_needed)];
		fill_chaotic_sorting_martrix(fsm_An, 0, 0, dimensions_needed, dimensions_needed * dimensions_needed, dimensions_needed, 1, fsm);
		cv::Mat img_copy = cv::Mat::zeros(image_to_resize.rows, image_to_resize.cols, CV_8UC3);
		cv::Mat blue_channel;
		cv::Mat green_channel;
		cv::Mat red_channel;
		cv::Mat channels[3];
		cv::Mat channels_orig[3];
		cv::split(img_copy, channels);
		cv::split(image_to_resize, channels_orig);
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
					channels[channel].at<uchar>(row_s1, col_s1) = channels_orig[channel].at<uchar>(row_fsm, col_fsm);
				}
			}
		}
		cv::merge(channels, 3, img_copy);
		FILE* fptr_write = fopen("keys.bin", "wb");
		fwrite(keys_array.data(), sizeof(double), 4, fptr_write);
		fclose(fptr_write);
		std::cout << std::endl;
		cv::imwrite("encrypted.png", img_copy);
		return true;
	}
};
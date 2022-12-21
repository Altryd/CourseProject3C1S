#pragma warning(disable : 4996)
#include <opencv2/opencv.hpp>
#include <openssl/sha.h>
#include <string_view>
#include <vector>
#include <array>
#include <fstream>
#include <algorithm>
#include <unordered_map>
#include <cmath>
#include <filesystem>
#include <chrono>
#include "ImageDecryptor.h"
#include "ImageEncryptor.h"

void print_info()
{
	std::cout << "Use --enc <path_to_file> <string_to_form_hash> to encrypt an image" << std::endl;
	std::cout << "Use --dec <path_to_file> --keys {path_to_keys} to decrypt an image" << std::endl;
}


int main(int argc, char** argv)
{
	setlocale(LC_ALL, "Ru");
	cv::Mat image;
	std::string string_to_hash;
	bool encryption = false;
	bool decryption = false;
	switch (argc)
{
	case 1:
	{
		print_info();
		return 0;
	}
	//TODO: argparse
	//TODO: openmp?
	//TODO: check chaotic encryption parallel algorithms
	case 4:
	{
		std::string argument(argv[1]);
		if (argument == "--enc")
		{
			encryption = true;
			std::string filepath(argv[2]);
			string_to_hash = argv[3];
			if (std::filesystem::exists(filepath))
			{
				image = cv::imread(filepath);
			}
			else
			{
				std::cout << "File don't exist, check path argument" << std::endl;
				return -1;
			}
		}
		else
		{
			print_info();
			return 0;
		}
		break;
	}
	case 5:
	{
		std::string dec_argument(argv[1]);
		std::string keys_argument(argv[3]);
		if (dec_argument == "--dec" && keys_argument == "--keys")
		{
			decryption = true;
			std::string filepath(argv[2]);
			std::string keys_path(argv[4]);
			if (std::filesystem::exists(filepath) && std::filesystem::exists(keys_path))
			{
				image = cv::imread(filepath);
			}
			else
			{
				std::cout << "File don't exist, check path argument" << std::endl;
				return -1;
			}
		}
		else
		{
			print_info();
			return 0;
		}
		break;
	}
	default:
	{
		print_info();
		return 0;
	}
}
	if (encryption)
	{
		/*auto begin = std::chrono::steady_clock::now();
		if (image.empty())
		{
			std::cout << "Image Not Found!!!" << std::endl;
			std::cin.get(); //wait for any key press
			return -1;
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
		std::string test_image = image_to_string(image);
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
		auto end = std::chrono::steady_clock::now();
		std::chrono::duration<double> diff = end - begin;
		std::cout << "The encryption took " << diff.count() << " seconds" << std::endl;
		std::cout << "The image has been successfully encrypted\nEncrypted image: encrypted.png\nKeys: keys.bin" << std::endl;*/
		ImageEncryptor encryptor(image, string_to_hash);
		try
		{
			auto begin = std::chrono::steady_clock::now();
			if (encryptor.encrypt_image())
			{
				auto end = std::chrono::steady_clock::now();
				std::chrono::duration<double> diff = end - begin;
				std::cout << "The encryption took " << diff.count() << " seconds" << std::endl;
				std::cout << "The image has been successfully encrypted\nEncrypted image: encrypted.png\nKeys: keys.bin" << std::endl;
			}
		}
		catch (const std::runtime_error& ex)
		{
			std::cout << ex.what() << std::endl;
			std::cin.get(); //wait for any key press
		}
		catch (...)
		{
			std::cout << "An unknown error has occured, try again" << std::endl;
		}
	}
	if (decryption)
	{
		ImageDecryptor decrypter(image, argv[4]);
		try
		{
			auto begin = std::chrono::steady_clock::now();
			if (decrypter.decrypt_image())
			{
				auto end = std::chrono::steady_clock::now();
				std::chrono::duration<double> diff = end - begin;
				std::cout << "The decryption took " << diff.count() << " seconds" << std::endl;
				std::cout << "The image has been successfully decrypted\nDecrypted image: decrypted.png" << std::endl;
			}
		}
		catch (const std::runtime_error& ex)
		{
			std::cout << ex.what() << std::endl;
			std::cin.get(); //wait for any key press
		}
		catch (...)
		{
			std::cout << "An unknown error has occured, try again" << std::endl;
		}
		/*
		auto begin = std::chrono::steady_clock::now();
		if (image.rows % 2 != 0 || image.cols != image.rows)
		{
			std::cout << "The image is corrupted" << std::endl;
			return -1;
		}
		int dimensions_needed = image.rows;
		std::string keys_path = argv[4];
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
		auto end = std::chrono::steady_clock::now();
		std::chrono::duration<double> diff = end - begin;
		std::cout << "The decryption took " << diff.count() << " seconds" << std::endl;
		std::cout << "The image has been successfully decrypted\nDecrypted image: decrypted.png" << std::endl;*/
	}

	return 0;
}
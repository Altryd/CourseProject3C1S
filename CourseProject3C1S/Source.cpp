#pragma warning(disable : 4996)
#include <opencv2/opencv.hpp>
#include <openssl/sha.h>
#include <string_view>
#include <vector>
#include <fstream>
#include <algorithm>
#include <unordered_map>

std::vector<int> sha256_ret_keys(const std::string& str)
{
	unsigned char hash[SHA256_DIGEST_LENGTH];
	SHA256_CTX sha256;
	SHA256_Init(&sha256);
	SHA256_Update(&sha256, str.c_str(), str.size());
	SHA256_Final(hash, &sha256);
	// const char* reinterpreted_const_unsigned = reinterpret_cast<const char*>(hash);
	// std::string_view string_view_test(reinterpreted_const_unsigned, sizeof(hash));
	// std::cout << string_view_test << std::endl;
	const short* reinterpreted_const_short = reinterpret_cast<const short*>(hash);
	const int* reinterpreted_const_int = reinterpret_cast<const int*>(hash);
	std::vector<int> keys;
	for (size_t i = 0; i < 6; i += 2)
	{
		keys.push_back(reinterpreted_const_short[i] / reinterpreted_const_short[i + 1]);
	}
	keys.push_back(reinterpreted_const_int[3] % 1000 + 1000);

	/*int k1 = reinterpreted_const_short[0] / reinterpreted_const_short[1];

	std::stringstream ss;
	for (int i = 0; i < SHA256_DIGEST_LENGTH; i++)
	{
		ss << std::hex << std::setw(2) << std::setfill('0') << (int)hash[i];
	} */
	return keys;
}


std::unordered_map<int, std::vector<int>> get_chaotic_sequences(const std::vector<int>& initial_keys, const int N)
{
	std::unordered_map<int, std::vector<int>> C_sequences;
	std::vector<int> values_past = initial_keys;
	values_past.pop_back();
	std::vector<int> values_present = values_past;
	for (size_t i = 0; i < static_cast<size_t>(initial_keys.back()); ++i)
	{
		values_present[0] = 35 * (values_past[1] - values_past[0]);  // x_new = 35*(y - x)
		values_present[1] = -7 * values_past[0] - values_past[0] * values_past[2] + 28 * values_past[1]; // y_new = -7*x - xz + 28y
		values_present[2] = values_past[0] * values_past[1] - 3 * values_past[2]; // z_new = x*y - 3*z
		values_past = values_present;
	}
	size_t iteration_count = static_cast<size_t>(N);
	iteration_count *= iteration_count;
	for (size_t i = 0; i < iteration_count; ++i)
	{
		values_present[0] = 35 * (values_past[1] - values_past[0]);  // x_new = 35*(y - x)
		values_present[1] = -7 * values_past[0] - values_past[0] * values_past[2] + 28 * values_past[1]; // y_new = -7*x - xz + 28y
		values_present[2] = values_past[0] * values_past[1] - 3 * values_past[2]; // z_new = x*y - 3*z
		values_past = values_present;
		C_sequences[1].push_back(values_present[0]);   //C1
		C_sequences[2].push_back(values_present[1]);   //C2
		C_sequences[3].push_back(values_present[2]);   //C3
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


std::string sha512(const std::string& str)
{
	unsigned char hash[SHA512_DIGEST_LENGTH];
	SHA512_CTX sha512;
	SHA512_Init(&sha512);
	SHA512_Update(&sha512, str.c_str(), str.size());
	SHA512_Final(hash, &sha512);
	std::stringstream ss;
	for (int i = 0; i < SHA512_DIGEST_LENGTH; i++)
	{
		ss << std::hex << std::setw(2) << std::setfill('0') << (int)hash[i];
	}
	return ss.str();
}

std::string image_to_string(const cv::Mat& img) {
	cv::Mat1b linear_img(img.reshape(1));
	return std::string(linear_img.begin(), linear_img.end());
}



void fill_chaotic_sorting_martrix(int* data, const int start_col, const int start_row,
	const int initial_cols, const int size, const int cols, const int start_val, const int fsm[4])
{
	if (size == 1)
	{
		//std::cout << start_row * initial_cols + start_col << "  writed: " << start_val << std::endl;
		data[start_row * initial_cols + start_col] = start_val;
		return;
	}
	fill_chaotic_sorting_martrix(data, start_col, start_row, initial_cols, size / 4, cols / 2, start_val + (fsm[0] - 1) * size / 4, fsm);
	fill_chaotic_sorting_martrix(data, start_col + cols / 2, start_row, initial_cols, size / 4, cols / 2, start_val + (fsm[1] - 1) * size / 4, fsm);
	fill_chaotic_sorting_martrix(data, start_col, start_row + cols / 2, initial_cols, size / 4, cols / 2, start_val + (fsm[2] - 1) * size / 4, fsm);
	fill_chaotic_sorting_martrix(data, start_col + cols / 2, start_row + cols / 2, initial_cols, size / 4, cols / 2, start_val + (fsm[3] - 1) * size / 4, fsm);
}



int main()
{
	setlocale(LC_ALL, "Ru");
	cv::Mat image = cv::imread("C:\\Users\\Altryd\\Downloads\\odd.png");  //TODO: change to relative path :)
	// image: Vec3b , а еще BGR здесь используется!
	// std::cout << cv::typeToString(image.type()) << std::endl;
	// Check for failure
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
	cv::Mat image_to_resize(dimensions_needed, dimensions_needed, image.type());
	cv::resize(image, image_to_resize, image_to_resize.size());
	for (int i = 0; i < image.rows; ++i)
	{
		for (int j = 0; j < image.cols; ++j)
		{
			image_to_resize.at<cv::Vec3b>(i, j) = cv::Vec3b(image.at<cv::Vec3b>(i, j));
		}
	}
	//   cv::imshow("Image Window Name here", image_to_resize);
	// Wait for any keystroke in the window
	//    cv::waitKey(0);
	std::string test_image = image_to_string(image);
	std::vector<int> sha256_res = sha256_ret_keys(test_image);
	std::cout << "KEYS: ";
	for (const auto it : sha256_res)
	{
		std::cout << it << std::endl;
	}
	auto sequences = get_chaotic_sequences(sha256_res, dimensions_needed);
	auto C1 = sequences[1];
	sequences[1].clear();
	auto S1 = sort_ranking(C1);
	C1.clear();

	int fsm[4] = { 4, 2, 1, 3 };
	int* fsm_An = new int[dimensions_needed * dimensions_needed];
	// std::cout << "dim needed:" << dimensions_needed;
	// std::cout << "dim needed square:" << dimensions_needed*dimensions_needed;
	fill_chaotic_sorting_martrix(fsm_An, 0, 0, dimensions_needed, dimensions_needed * dimensions_needed, dimensions_needed, 1, fsm);

	std::ofstream myfile("example.txt");
	if (myfile.is_open())
	{
		for (size_t i = 0; i < dimensions_needed; ++i)
		{
			for (size_t j = 0; j < dimensions_needed; ++j)
			{
				myfile << fsm_An[i * dimensions_needed + j] << "  ";
			}
			myfile << std::endl;
		}
	}

	cv::imshow("orig", image_to_resize);
	cv::waitKey(0);
	auto img_copy = image_to_resize.clone();  //CRUCIAL !!
	cv::imshow("orig2_withclone", image_to_resize);
	cv::waitKey(0);
	for (size_t i = 0; i < dimensions_needed; ++i)
	{
		for (size_t j = 0; j < dimensions_needed; ++j)
		{
			auto row_s1 = S1[i * dimensions_needed + j] / dimensions_needed;
			auto col_s1 = S1[i * dimensions_needed + j] % dimensions_needed;
			auto row_fsm = (fsm_An[i * dimensions_needed + j] - 1) / dimensions_needed;
			auto col_fsm = (fsm_An[i * dimensions_needed + j] - 1) % dimensions_needed;
			img_copy.at<cv::Vec3b>(row_s1, col_s1) = image_to_resize.at<cv::Vec3b>(row_fsm, col_fsm);
		}
	}
	cv::imshow("half-encrypted", img_copy);
	auto decrypted = img_copy.clone();
	cv::waitKey(0);
	for (size_t i = 0; i < dimensions_needed; ++i)
	{
		for (size_t j = 0; j < dimensions_needed; ++j)
		{
			auto row_s1 = S1[i * dimensions_needed + j] / dimensions_needed;
			auto col_s1 = S1[i * dimensions_needed + j] % dimensions_needed;
			auto row_fsm = (fsm_An[i * dimensions_needed + j] - 1) / dimensions_needed;
			auto col_fsm = (fsm_An[i * dimensions_needed + j] - 1) % dimensions_needed;
			decrypted.at<cv::Vec3b>(row_fsm, col_fsm) = img_copy.at<cv::Vec3b>(row_s1, col_s1);
		}
	}
	cv::imshow("decrypted?", decrypted);
	cv::waitKey(0);



	auto C2 = sequences[2];
	std::vector<int> T(dimensions_needed * dimensions_needed);
	for (size_t i = 0; i < T.size(); ++i)
	{
		T[i] = (int)ceil(C2[i] * pow(10, 5)) % 256;  // TODO: округление вверх
	}
	delete[] fsm_An;
	// auto sequences = get_chaotic_sequences(sha256_res, dimensions_needed);


	/*TODO: #5 task from habr:
	auto C3 = sequences[3];
	sequences[3].clear();
	auto S3 = sort_ranking(C3);
	C3.clear();
	auto P1 = linear_img;
	P1[0] = T[0] ^ linear_img.at<cv::Vec3b>(0);
	for (size_t i = 1; i < dimensions_needed * dimensions_needed; ++i)
	{
		P1[i] = P1[i - 1] ^ T[i] ^ linear_img.at<cv::Vec3b>(i);
	}*/


	return 0;
}
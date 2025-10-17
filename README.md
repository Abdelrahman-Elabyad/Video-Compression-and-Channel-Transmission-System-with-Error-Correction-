# Video Transmission with Convolutional and Punctured Coding

This MATLAB project simulates **video transmission over a noisy channel** using various **convolutional and punctured coding schemes**.  
The experiment evaluates **Bit Error Rate (BER)** and **Throughput performance** under different channel error probabilities.

---

## 🎯 Project Overview
The system transmits the `highway.avi` video through a **Binary Symmetric Channel (BSC)** with different bit error probabilities (`p`).  
Each frame of the video is:
1. Converted into binary data streams.  
2. Encoded using multiple **convolutional** and **punctured codes**.  
3. Transmitted through a simulated noisy channel.  
4. Decoded using **Viterbi decoding**.  
5. Reconstructed into video frames for visual comparison.

---

## ⚙️ Implemented Coding Schemes
- **Uncoded Transmission** (baseline)  
- **Convolutional Codes**  
  - Rate 1/2  
  - Rate 1/3  
  - Rate 1/4  
  - Rate 2/3  
- **Punctured Codes**  
  - Rate 8/9  
  - Rate 4/5  
  - Rate 2/3  

Each scheme is compared across a range of error probabilities from **0.001 to 0.2**.

---

## 📈 Results and Analysis
- The program computes and plots:
  - **Bit Error Rate (BER)** vs. Channel Error Probability  
  - **Throughput (Effective Rate)** vs. Channel Error Probability  
- The generated plots allow for clear performance comparison among all coding techniques.
- Output videos are also generated for extreme cases (`p = 0.001` and `p = 0.1`) to visualize the degradation in quality.

---

## 🧩 Key MATLAB Functions Used
- `VideoReader`, `VideoWriter` – for video input/output  
- `poly2trellis`, `convenc`, `vitdec` – for convolutional coding and decoding  
- `bsc` – for binary symmetric channel simulation  
- `bi2de`, `de2bi`, `reshape` – for binary-pixel conversions  

---

## 📦 Output Files
For each tested probability:
- BER and Throughput plots are saved.
- Reconstructed videos are saved for:
  - `no_coding`
  - `conv_1_2`, `conv_1_3`, `conv_1_4`, `conv_2_3`
  - `punc_8_9`, `punc_4_5`, `punc_2_3`

Example output file:
no_coding_p0.001000.avi
conv_1_2_p0.001000.avi

---

## 🧠 Purpose
This project demonstrates the **effectiveness of error-control coding** in multimedia transmission over noisy communication channels.  
It shows how **coding rate, puncturing, and noise level** impact the balance between **error resilience** and **throughput efficiency**.

---

## 👨‍💻 Author
Developed for educational and research purposes in **Digital Communication Systems** and **Channel Coding** coursework.

---

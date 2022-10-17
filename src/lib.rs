use std::f64::consts::PI;

// generate complex phasor for FFT computation
pub fn generate_wn(sample_len: usize) -> Vec<(f64, f64)> {
    let sample_len_f: f64 = sample_len as f64;
    let w: f64 = 2.0 * PI / sample_len_f;
    let mut wn: Vec<(f64, f64)> = vec![(0.0, 0.0); sample_len / 2];
    let mut n: f64;
    let mut phi: f64;
    
    for index in 0..(sample_len / 2) {
        n = index as f64;        
        phi = w * n;
        wn[index].0 = phi.cos();        
        wn[index].1 = -phi.sin();                    
    }   
     
    wn
}

// main FFT algorithm, real input signals only
//
// splits DFT Y(k) into two DFTs of half the size, G(k) and H(k), 
// where G(k) takes the even samples and H(k) takes the odd samples of input to Y(k).
// the outputs of original DFT are reconstructed given Y(k) = G(k) + H(k) * Wn(k),
// where Wn(k) is a complex phasor with frequency of the sample period.
// 
// if log2(SAMPLE_LEN) is an integer then DFTs G(k) and H(k) can be further split into DFTs of half size.
// this process repeats until the number of DFTs equals the sample length, each with 1 input/ouput.
// for a 1 input DFT, Y(0) = x(0). therefore, no DFTs need to be computed, only reconstructed with Wn.
pub fn fft(input: Vec<f64>, wn: Vec<(f64, f64)>, sample_len: usize) -> Vec<(f64, f64)> { 
    let sample_len_int = sample_len as i64;
    let sample_len_f: f64 = sample_len as f64;
    let num_bits = sample_len_f.log2().round() as i64;
    let mut g_index: usize;
    let mut h_index: usize;
    let mut wn_index: usize;
    let mut dft_size: i64 = 2;
    let mut dft_num: i64;
    let mut x_real = reverse_bin_index(input, num_bits, sample_len);
    let mut x_imag: Vec<f64> = vec![0.0; sample_len];
    let mut y_real: Vec<f64> = vec![0.0; sample_len];
    let mut y_imag: Vec<f64> = vec![0.0; sample_len];
    let mut h_wn_real: f64;
    let mut h_wn_imag: f64;
    
    while dft_size <= sample_len_int {
        dft_num = sample_len_int / dft_size;

        // iterate through DFTs
        for dft_id in 0..dft_num {

            // iterate through DFT indices
            for k in 0..(dft_size / 2) {
                g_index = (k + dft_id * dft_size) as usize;                // G(k) index
                h_index = (k + dft_id * dft_size + dft_size / 2) as usize; // H(k) index
                wn_index = (k * dft_num) as usize;                         // current iteration takes even Wn(k) indices of next iteration
                                
                // calculate H(k) * Wn(k)
                h_wn_real = x_real[h_index] * wn[wn_index].0 - x_imag[h_index] * wn[wn_index].1;
                h_wn_imag = x_real[h_index] * wn[wn_index].1 + x_imag[h_index] * wn[wn_index].0;

                // calculate current DFT output Y(k) = G(k) + H(k) * Wn(k)
                y_real[g_index] = x_real[g_index] + h_wn_real;
                y_imag[g_index] = x_imag[g_index] + h_wn_imag;
                y_real[h_index] = x_real[g_index] - h_wn_real;  // due to symmetry of complex phasor Wn(k) and periodicity of G(k) and H(k)
                y_imag[h_index] = x_imag[g_index] - h_wn_imag;  // when k > dft_size / 2,
                                                                // Y(k) = G(k - dft_size / 2) + H(k - dft_size / 2) * -Wn(k - dft_size / 2)                
            }
        }
        
        // double size of DFT and set input to output values for next iteration
        dft_size *= 2;
        for n in 0..sample_len {
            x_real[n] = y_real[n];
            x_imag[n] = y_imag[n];
        }
    }
    
    // convert output to polar coordinates
    cart_to_polar(y_real, y_imag, sample_len)
}

// convert list of complex numbers from cartesian to polar coordinates
fn cart_to_polar(input_real: Vec<f64>, input_imag: Vec<f64>, sample_len: usize) -> Vec<(f64, f64)> {
    let mut polar_list: Vec<(f64, f64)> = vec![(0.0, 0.0); sample_len];
    
    for n in 0..sample_len {
        polar_list[n].0 = complex_mag(input_real[n], input_imag[n]);
        polar_list[n].1 = complex_phase(input_real[n], input_imag[n]);
    }

    polar_list   
}

// calculate magnitude of complex number
fn complex_mag(a: f64, b: f64) -> f64 {
    (a.powf(2.0) + b.powf(2.0)).sqrt()
}

// calculate phase of complex number
fn complex_phase(a: f64, b: f64) -> f64 {
    b.atan2(a)   
}

// reverse bits of binary number
fn reverse_bits(mut input: i64, num_bits: i64) -> i64 {
    let mut output: i64 = 0;

    for _ in 0..num_bits {
        output <<= 1;
        if input & 1 == 1 {
            output ^= 1;
        }       
        input >>= 1;
    }
    
    output
}

// reorder list so bits of binary index are reversed
fn reverse_bin_index(input: Vec<f64>, num_bits: i64, sample_len: usize) -> Vec<f64> {
    let mut output: Vec<f64> = vec![0.0; sample_len];
    let mut n: i64;
    let mut rev_index: usize;
        
    for index in 0..sample_len {
        n = index as i64;
        rev_index = reverse_bits(n, num_bits) as usize;
        output[index] = input[rev_index];
    }
    
    output
}

// check that log2(sample length) is an integer
fn check_sample_len(sample_len: f64) -> bool {
    let y: f64 = sample_len.log2();

    if y % 1.0 < 1e-10 {
        true
    }
    else {
        println!("error: log2(sample length) not an integer");
        false
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::Instant;
    
    // generate sinusoid with harmonic of sample period for creating test signals
    // fundamental frequency is harmonic = 1
    fn generate_sinusoid(harmonic: i64, mag: f64, phase: f64, sample_len: usize) -> Vec<f64> {
        let sample_len_f = sample_len as f64;
        let w: f64 = 2.0 * PI / sample_len_f;
        let harmonic_f = harmonic as f64;
        let mut output: Vec<f64> = vec![0.0; sample_len];
        let mut n: f64;
        let mut phi: f64;
    
        for index in 0..sample_len {
            n = index as f64;        
            phi = w * harmonic_f * n + phase;
            output[index] = mag * phi.cos();        
        }   
    
        output
    }

    // generate a test signal for FFT input
    fn test_signal(sample_len: usize) -> Vec<f64> {
        let sig1: Vec<f64> = generate_sinusoid(1, 1.0, 0.0, sample_len);
        let sig2: Vec<f64> = generate_sinusoid(2, 0.5, PI, sample_len);
        let mut output: Vec<f64> = vec![0.0; sample_len];
    
        for n in 0..sample_len {
            output[n] = sig1[n] + sig2[n];
        }
    
        output
    }

    // inverse DFT for testing FFT output
    fn inv_dft(input: Vec<(f64, f64)>, sample_len: usize) -> Vec<f64> {
        let sample_len_f = sample_len as f64;
        let w: f64 = 2.0 * PI / sample_len_f;
        let mut output: Vec<f64> = vec![0.0; sample_len];
        let mut n: f64;
        let mut k: f64;
        let mut phi: f64;
    
        for n_index in 0..sample_len {
            n = n_index as f64;        
        
            for k_index in 0..sample_len {
                k = k_index as f64;

                phi = w * n * k + input[k_index].1;
                output[n_index] += input[k_index].0 * phi.cos();
            }
              
            output[n_index] /= sample_len_f;        
        }   
    
        output
    }

    // DFT for testing inverse DFT
    fn dft(input: Vec<f64>, sample_len: usize) -> Vec<(f64, f64)> {
        let sample_len_f = sample_len as f64;
        let w: f64 = 2.0 * PI / sample_len_f;
        let mut y_real: Vec<f64> = vec![0.0; sample_len];
        let mut y_imag: Vec<f64> = vec![0.0; sample_len];
        let mut n: f64;
        let mut k: f64;
        let mut phi: f64;
    
        for k_index in 0..sample_len {
            k = k_index as f64;        
        
            for n_index in 0..sample_len {
                n = n_index as f64;

                phi = w * n * k;
                y_real[k_index] += input[n_index] * phi.cos();
                y_imag[k_index] += input[n_index] * -phi.sin();
            }             
        }   
    
        cart_to_polar(y_real, y_imag, sample_len)
    }

    // time FFT calculation time of increasingly large sample lengths
    #[test]
    fn speed_test(){
        let mut num_bits: u32 = 1;
        let mut sample_len: usize;
        
        while num_bits <= 16 {
            sample_len = 2usize.pow(num_bits);
            println!("\nsample length: {}", sample_len);
            println!("number of bits: {}", num_bits);

            let wn: Vec<(f64, f64)> = generate_wn(sample_len);   
            let dft_input: Vec<f64> = test_signal(sample_len);   
            let fft_input: Vec<f64> = dft_input.clone();
        
            // compare time of FFT algorithm with DFT
            //let dft_timer = Instant::now();         
            //let _output_dft: Vec<(f64, f64)> = dft(dft_input, sample_len);
            //println!("DFT Elapsed time: {:.2?}", dft_timer.elapsed());

            let fft_timer = Instant::now();         
            let _output_fft: Vec<(f64, f64)> = fft(fft_input, wn, sample_len);
            println!("FFT Elapsed time: {:.2?}", fft_timer.elapsed());
            
            num_bits += 1;
        }
    }

    // test FFT by recombining FFT co-efficients into original signal
    #[test]
    fn fft_test() {
        let mut e: f64;
        let sample_len: usize = 32;
        let wn: Vec<(f64, f64)> = generate_wn(sample_len);   
        println!("\nWn: {:?}", wn);
        let input: Vec<f64> = test_signal(sample_len);         
        let input_clone: Vec<f64> = input.clone();
        println!("\ninput: {:?}", input);
        let output: Vec<(f64, f64)> = fft(input, wn, sample_len);   
        println!("\noutput: {:?}", output);   
        let recombination: Vec<f64> = inv_dft(output, sample_len);
        println!("\nrecombination: {:?}", recombination);

        for n in 0..sample_len {
            e = recombination[n] - input_clone[n];
            assert!(e < 1e-10);
        }
    }

    // test inverse DFT by performing DFT on a signal and recombining
    #[test]
    fn inv_dft_test() {  
        let mut e: f64;
        let sample_len: usize = 32;
        let input: Vec<f64> = test_signal(sample_len);
        let input_clone: Vec<f64> = input.clone();
        println!("\ninput: {:?}", input);
        let output: Vec<(f64, f64)> = dft(input, sample_len);    
        println!("\noutput: {:?}", output);
        let recombination: Vec<f64> = inv_dft(output, sample_len); 
        println!("\nrecombination: {:?}", recombination);
    
        for n in 0..sample_len {
            e = recombination[n] - input_clone[n];
            assert!(e < 1e-10);
        }
    }

    // tests bit reversal for all sample lengths up to 64-bit
    #[test]
    fn reverse_bits_test() {
        let mut input: i64;
        let mut output: i64;
        let mut input_bit: i64;
        let mut output_bit: i64;
        let mut sample_len: usize = 2;
        let mut sample_len_f = sample_len as f64;
        let mut num_bits = sample_len_f.log2().round() as i64;    

        while num_bits <= 62 {
            println!("\nsample length: {}", sample_len);
            println!("number of bits: {}", num_bits);

            input = 1;
            
            for index in 0..num_bits {
                output = reverse_bits(input, num_bits);       

                //println!("\ninput: {}", input);
                //println!("output: {}", output);   

                for bit in 0..num_bits {
                    input_bit = (input >> (num_bits - bit - 1)) & 1;
                    output_bit = (output >> bit) & 1;
                    assert_eq!(input_bit, output_bit);
                }
                
                input *= 2;
            }  

            sample_len *= 2;
            sample_len_f = sample_len as f64;
            num_bits = sample_len_f.log2().round() as i64;
        }
    }
}
use std::f64::consts::PI;
use std::time::Instant;

pub const SAMPLE_LEN: usize = 16;
pub const SAMPLE_LEN_F: f64 = SAMPLE_LEN as f64;
pub const W: f64 = 2.0 * PI / SAMPLE_LEN_F;

// generate complex phasor for FFT computation
fn generate_wn() -> [(f64, f64); SAMPLE_LEN / 2] {
    let mut wn: [(f64, f64); SAMPLE_LEN / 2] = [(0.0, 0.0); SAMPLE_LEN / 2];
    let mut n: f64;
    let mut phi: f64;
    
    for index in 0..(SAMPLE_LEN / 2) {
        n = index as f64;        
        phi = W * n;
        wn[index].0 = phi.cos();        
        wn[index].1 = -phi.sin();                    
    }   
     
    wn
}

// main FFT algorithm 
// real input signals only, for complex inputs add imaginary input as parameter and set x_imag.
//
// splits DFT Y(k) into two DFTs of half the size, G(k) and H(k), 
// where G(k) takes the even samples and H(k) takes the odd samples of input to Y(k).
// the outputs of original DFT are reconstructed given Y(k) = G(k) + H(k) * Wn(k),
// where Wn(k) is a complex phasor with frequency of the sample period.
// 
// if log2(SAMPLE_LEN) is an integer then DFTs G(k) and H(k) can be further split into DFTs of half size.
// this process repeats until the number of DFTs equals the sample length, each with 1 input/ouput.
// for a 1 input DFT, Y(0) = x(0). therefore, no DFTs need to be computed, only reconstructed with Wn.
fn fft(input: [f64; SAMPLE_LEN], wn: [(f64, f64); SAMPLE_LEN / 2]) -> [(f64, f64); SAMPLE_LEN] { 
    let mut g_index: usize;
    let mut h_index: usize;
    let mut wn_index: usize;
    let sample_len_int = SAMPLE_LEN as i64;
    let mut dft_size: i64 = 2;
    let mut dft_num: i64;
    let mut x_real = input;
    let mut x_imag: [f64; SAMPLE_LEN] = [0.0; SAMPLE_LEN];
    let mut y_real: [f64; SAMPLE_LEN] = [0.0; SAMPLE_LEN];
    let mut y_imag: [f64; SAMPLE_LEN] = [0.0; SAMPLE_LEN];
    let mut h_wn_real: f64;
    let mut h_wn_imag: f64;
    
    while dft_size <= sample_len_int {
        dft_num = sample_len_int / dft_size;
        //println!("\nDFT size: {}", dft_size);

        // iterate through DFTs
        for dft_id in 0..dft_num {
            //println!("\nDFT id: {}", dft_id);

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
                //println!("\nDFT index: {}", k);
                //println!("g_index: {}", g_index);
                //println!("h_index: {}", h_index);
                //println!("wn_index: {}", wn_index);
                //println!("\nx_real[h_index]: {}", x_real[h_index]);
                //println!("x_imag[h_index]: {}", x_imag[h_index]);
                //println!("wn[wn_index]: {} {}", wn[wn_index].0, wn[wn_index].1);
                //println!("\ny_real: {:?}", y_real);
                //println!("\ny_imag: {:?}", y_imag);
            }
        }
        
        // double size of DFT and set input to output values for next iteration
        dft_size *= 2;
        for n in 0..SAMPLE_LEN {
            x_real[n] = y_real[n];
            x_imag[n] = y_imag[n];
        }
    }
    
    // convert output to polar coordinates
    cart_to_polar(y_real, y_imag)
}

// convert list of complex numbers from cartesian to polar coordinates
fn cart_to_polar(input_real: [f64; SAMPLE_LEN], input_imag: [f64; SAMPLE_LEN]) -> [(f64, f64); SAMPLE_LEN] {
    let mut polar_list: [(f64, f64); SAMPLE_LEN] = [(0.0, 0.0); SAMPLE_LEN];
    
    for n in 0..SAMPLE_LEN {
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
fn reverse_bin_index(input: [f64; SAMPLE_LEN], num_bits: i64) -> [f64; SAMPLE_LEN] {
    let mut output: [f64; SAMPLE_LEN] = [0.0; SAMPLE_LEN];
    let mut n: i64;
    let mut rev_index: usize;
        
    for index in 0..SAMPLE_LEN {
        n = index as i64;
        rev_index = reverse_bits(n, num_bits) as usize;
        output[index] = input[rev_index];
    }
    
    output
}

// check that log2(sample length) is an integer
fn check_sample_len() -> bool {
    let y: f64 = SAMPLE_LEN_F.log2();

    if y % 1.0 < 1e-10 {
        true
    }
    else {
        println!("error: log2(sample length) not an integer");
        false
    }
}



// generate sinusoid with harmonic of sample period for creating test signals
// fundamental frequency is harmonic = 1
fn generate_sinusoid(harmonic: i64, mag: f64, phase: f64) -> [f64; SAMPLE_LEN] {
    let harmonic_f: f64 = harmonic as f64;
    let mut output: [f64; SAMPLE_LEN] = [0.0; SAMPLE_LEN];
    let mut n: f64;
    let mut phi: f64;
    
    for index in 0..SAMPLE_LEN {
        n = index as f64;        
        phi = W * harmonic_f * n + phase;
        output[index] = mag * phi.cos();        
    }   
    
    output
}

// generate a test signal for FFT input
fn test_signal() -> [f64; SAMPLE_LEN] {
    let sig1: [f64; SAMPLE_LEN] = generate_sinusoid(1, 1.0, 0.0);
    let sig2: [f64; SAMPLE_LEN] = generate_sinusoid(2, 0.5, PI);
    let mut output: [f64; SAMPLE_LEN] = [0.0; SAMPLE_LEN];
    
    for n in 0..SAMPLE_LEN {
        output[n] = sig1[n] + sig2[n];
    }
    
    output
}

// inverse DFT for testing FFT output
fn inv_dft(input: [(f64, f64); SAMPLE_LEN]) -> [f64; SAMPLE_LEN] {
    let mut output: [f64; SAMPLE_LEN] = [0.0; SAMPLE_LEN];
    let mut n: f64;
    let mut k: f64;
    let mut phi: f64;
    
    for n_index in 0..SAMPLE_LEN {
        n = n_index as f64;        
        
        for k_index in 0..SAMPLE_LEN {
            k = k_index as f64;

            phi = W * n * k + input[k_index].1;
            output[n_index] += input[k_index].0 * phi.cos();
        }
              
        output[n_index] /= SAMPLE_LEN_F;        
    }   
    
    output
}

// DFT for testing inverse DFT
fn dft(input: [f64; SAMPLE_LEN]) -> [(f64, f64); SAMPLE_LEN] {
    let mut y_real: [f64; SAMPLE_LEN] = [0.0; SAMPLE_LEN];
    let mut y_imag: [f64; SAMPLE_LEN] = [0.0; SAMPLE_LEN];
    let mut n: f64;
    let mut k: f64;
    let mut phi: f64;
    
    for k_index in 0..SAMPLE_LEN {
        k = k_index as f64;        
        
        for n_index in 0..SAMPLE_LEN {
            n = n_index as f64;

            phi = W * n * k;
            y_real[k_index] += input[n_index] * phi.cos();
            y_imag[k_index] += input[n_index] * -phi.sin();
        }             
    }   
    
    cart_to_polar(y_real, y_imag)
}

// public function for executing FFT lib. input not parameter yet.
pub fn execute_fft(input: [f64; SAMPLE_LEN]) -> [(f64, f64); SAMPLE_LEN] {
    if check_sample_len() {
        let num_bits = SAMPLE_LEN_F.log2() as i64;
        let wn: [(f64, f64); SAMPLE_LEN / 2] = generate_wn();   
        //let mut input: [f64; SAMPLE_LEN] = test_signal();   

        let x: [f64; SAMPLE_LEN] = reverse_bin_index(input, num_bits);
        fft(x, wn)
    }
    else {
        [(0.0, 0.0); SAMPLE_LEN]
    }   
}

#[test]
fn speed_test(){
    if check_sample_len() {
        let num_bits = SAMPLE_LEN_F.log2() as i64;
        let wn: [(f64, f64); SAMPLE_LEN / 2] = generate_wn();   
        let mut input: [f64; SAMPLE_LEN] = test_signal();   

        // compare time of FFT algorithm with DFT
        let time = Instant::now();         
        let _output_dft: [(f64, f64); SAMPLE_LEN] = dft(input);
        println!("DFT Elapsed time: {:.2?}", time.elapsed());

        let time = Instant::now();         
        input = reverse_bin_index(input, num_bits);
        let _output: [(f64, f64); SAMPLE_LEN] = fft(input, wn);
        println!("FFT Elapsed time: {:.2?}", time.elapsed());
        //println!("\noutput: {:?}", _output);
    }
}

#[test]
fn fft_test() {
    let mut e: f64;

    if check_sample_len() {
        let num_bits = SAMPLE_LEN_F.log2() as i64;
        let wn: [(f64, f64); SAMPLE_LEN / 2] = generate_wn();   
        let input: [f64; SAMPLE_LEN] = test_signal();         
        let input_rev: [f64; SAMPLE_LEN] = reverse_bin_index(input, num_bits);
        let output: [(f64, f64); SAMPLE_LEN] = fft(input_rev, wn);   
        let recombination: [f64; SAMPLE_LEN] = inv_dft(output);

        println!("\nWn: {:?}", wn);
        println!("\ninput: {:?}", input);
        println!("\ninput (re-ordered): {:?}", input_rev);
        println!("\noutput: {:?}", output);   
        println!("\nrecombination: {:?}", recombination);
    
        for n in 0..SAMPLE_LEN {
            e = recombination[n] - input[n];
            assert!(e < 1e-10);
        }
    }
}

#[test]
fn inv_dft_test() {  
    let mut e: f64;
    let input: [f64; SAMPLE_LEN] = test_signal();
    let output: [(f64, f64); SAMPLE_LEN] = dft(input);    
    let recombination: [f64; SAMPLE_LEN] = inv_dft(output); 

    println!("\ninput: {:?}", input);
    println!("\noutput: {:?}", output);
    println!("\nrecombination: {:?}", recombination);
    
    for n in 0..SAMPLE_LEN {
        e = recombination[n] - input[n];
        assert!(e < 1e-10);
    }
}

#[test]
fn reverse_bits_test() {
    let mut input: i64;
    let mut output: i64;
    let mut input_bit: i64;
    let mut output_bit: i64;
    let num_bits = SAMPLE_LEN_F.log2() as i64;
    
    println!("\nsample length: {}", SAMPLE_LEN);
    println!("number of bits: {}", num_bits);

    for index in 0..SAMPLE_LEN {
        input = index as i64;
        output = reverse_bits(input, num_bits);       

        println!("\ninput: {}", input);
        println!("output: {}", output);   

        for bit in 0..num_bits {
            input_bit = (input >> (num_bits - bit - 1)) & 1;
            output_bit = (output >> bit) & 1;
            assert_eq!(input_bit, output_bit);
        }
    }  
}

use fft::{execute_fft, SAMPLE_LEN};
use std::f64::consts::PI;  
pub const SAMPLE_LEN_F: f64 = SAMPLE_LEN as f64;                  
pub const W: f64 = 2.0 * PI / SAMPLE_LEN_F;

// generate sinusoid with harmonic of sample period for creating t
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

fn main() {
    let input: [f64; SAMPLE_LEN] = test_signal();
    let output: [(f64, f64); SAMPLE_LEN] = execute_fft(input);
    println!("main output: {:?}", output);
}


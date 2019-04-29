data { 
  int num_data; 
  int num_basis; 
  vector[num_data] Y; 
  matrix[num_data, num_basis] B; 
} 
 
parameters { 
  vector[num_basis] w; 
  vector<lower=0>[num_basis] w_sigma; 
} 
 
transformed parameters { 
  vector[num_data] mu;
  vector[num_data] sigma;
  mu = to_vector(B*w); 
  sigma = to_vector(B*w_sigma); 
} 
 
model { 
  w ~ normal(0, 1); 
  w_sigma ~ exponential(1); 
  Y ~ normal(mu, sigma); 
}


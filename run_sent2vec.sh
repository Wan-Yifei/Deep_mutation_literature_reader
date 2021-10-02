export MODEL_LOCATION=/home/ec2-user/BioSentVec_PubMed_MIMICIII-bigram_d700.bin
export RAW_CONTENT_PATH=/home/ec2-user/Deep_mutation_literature_reader/data/cleaned_abstract.pickle
export HMAC_KEY=wanyifei

python3 sent2vec_pipeline.py embedded_abstract.pickle

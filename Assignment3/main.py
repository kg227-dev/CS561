from part1 import *


if __name__ == '__main__':
    # Parse FASTA file input (incl. one-hot-encodings for sequences)

    sim1_data = parse_data("/Users/sydneyballard/Desktop/Desktop - Sydney’s MacBook Pro/CS 561/cs561 repository COLLABORATIVE W KUSH/CS561/Assignment3/sim1/test.fasta")
    sim2_data = parse_data("/Users/sydneyballard/Desktop/Desktop - Sydney’s MacBook Pro/CS 561/cs561 repository COLLABORATIVE W KUSH/CS561/Assignment3/sim2/test.fasta")
    sim6_data = parse_data("/Users/sydneyballard/Desktop/Desktop - Sydney’s MacBook Pro/CS 561/cs561 repository COLLABORATIVE W KUSH/CS561/Assignment3/sim6/test.fasta")
    sim7_data = parse_data("/Users/sydneyballard/Desktop/Desktop - Sydney’s MacBook Pro/CS 561/cs561 repository COLLABORATIVE W KUSH/CS561/Assignment3/sim7/test.fasta")

    sim_datasets = [sim1_data,
                    sim2_data,
                    sim6_data,
                    sim7_data]
    
    for dataset in sim6_data:
        # Make model

        # Train model
            model.fit(X_train, Y_train, verbose=1, validation_data=(X_valid, Y_valid),batch_size=128, epochs=100, callbacks=[EarlyStopping(patience=10, monitor="val_loss", restore_best_weights=True), History()])


        # Test model accuracy
            pred = model.predict(X_test, batch_size=128)
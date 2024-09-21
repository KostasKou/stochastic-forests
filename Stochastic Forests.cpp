	#include <iostream>
	#include <vector>
	#include <map>
	#include <algorithm>
	#include <unordered_map>
	#include <random>
	#include <iomanip>
	#include <fstream>

	using namespace std;
	struct DataPoint {
		vector<int> categoricalAttributes;
		vector<double> continuousAttributes;
		int target;
	};
	void split(const string& s, char c,vector<string>& v) {
	   int i = 0;
	   int j = s.find(c);
	   while (j >= 0) {
		  v.push_back(s.substr(i, j-i));
		  i = ++j;
		  j = s.find(c, j);
		  if (j < 0) {
			 v.push_back(s.substr(i, s.length()));
		  }
	   }
	}
	void loadCSV(istream& in,vector<DataPoint>& dataset) {
	   vector<string> types; //float,float,int,int,...
	   vector<string> p ; //vector 1 parathrhshs me teleutaio stoixeio to target.
	   string tmp;
	   int catAt;
	   double contAt;
	   int counter = 0;
	   while (!in.eof()) {
		  counter++;
		  DataPoint data;


		  if(counter == 1){
			getline(in, tmp, '\n');
			split(tmp, ',', p);
			for(int i=0; i<p.size(); i++){
				types.push_back(p[i]);
			}
			tmp.clear();
			p.clear();
		  }
		  else{
			getline(in, tmp, '\n');
			split(tmp, ',', p);
			for(int i=0; i<p.size()-1; i++){
				if(types[i]=="int"){
					catAt = stoi(p[i]);
					data.categoricalAttributes.push_back(catAt);

				}
				else if(types[i]=="float"){
					contAt = stod(p[i]);
					data.continuousAttributes.push_back(contAt);
					//cout<<"stod(pi)"<<stod(p[i])<<endl;
				}
			}
			data.target = stoi(p[p.size()-1]);
			dataset.push_back(data);
			tmp.clear();
			p.clear();

		  }
	   }
	}


	struct Node {
		bool leafNode;
		int leafValue;
		double splitValue;
		int attributeIndex;
		
		vector<Node*> children;

		Node() : leafNode(false), leafValue(0), attributeIndex(-1), splitValue(0.0) {}

		~Node() {
			for (auto& child : children) {
				delete child;
			}
		}
	};

	double calculateEntropy(const vector<DataPoint>& dataset) {
		map<int, int> classCounts;
		int numExamples = dataset.size();

		for (const auto& data : dataset) {
			classCounts[data.target]++;
		}

		double entropy = 0.0;

		for (const auto& entry : classCounts) {
			double probability = static_cast<double>(entry.second) / numExamples;
			entropy -= probability * log2(probability);
		}

		return entropy;
	}
	double getBestSplitValue(const vector<DataPoint>& dataset, int attributeIndex);
	// Calculate information gain for a categorical attribute split
	double calculateGain(const vector<DataPoint>& dataset, int attributeIndex, bool isCategorical,double splitValue) {
		double parentEntropy = calculateEntropy(dataset);
		int numExamples = dataset.size();
		double totalGain = 0.0;

		if (isCategorical) {
			map<int, vector<DataPoint>> attributeValues;

			for (const auto& data : dataset) {
				attributeValues[data.categoricalAttributes[attributeIndex]].push_back(data);
			}

			for (const auto& entry : attributeValues) {
				double subsetEntropy = calculateEntropy(entry.second);
				double subsetProbability = static_cast<double>(entry.second.size()) / numExamples;
				totalGain += subsetProbability * subsetEntropy;
			}
		} else { // Continuous attribute
			// Get the best split value for the continuous attribute
			splitValue = getBestSplitValue(dataset, attributeIndex);

			vector<DataPoint> leftSubset;
			vector<DataPoint> rightSubset;

			for (const auto& data : dataset) {
				if (data.continuousAttributes[attributeIndex] <= splitValue) {
					leftSubset.push_back(data);
				} else {
					rightSubset.push_back(data);
				}
			}

			double leftSubsetEntropy = calculateEntropy(leftSubset);
			double rightSubsetEntropy = calculateEntropy(rightSubset);
			double leftSubsetProbability = static_cast<double>(leftSubset.size()) / numExamples;
			double rightSubsetProbability = static_cast<double>(rightSubset.size()) / numExamples;

			totalGain = (leftSubsetProbability * leftSubsetEntropy) + (rightSubsetProbability * rightSubsetEntropy);
		}

		double gain = parentEntropy - totalGain;
		return gain;
	}


	// Get the best split value for a continuous attribute
	double calculateEntropy2(const std::unordered_map<int, int>& classCounts, int totalSize) {
		double entropy = 0.0;
		for (const auto& pair : classCounts) {
			double probability = static_cast<double>(pair.second) / totalSize;
			entropy -= probability * log2(probability);
		}
		return entropy;
	}

	double getBestSplitValue(const std::vector<DataPoint>& dataset, int attributeIndex) {
		// Initialize variables
		double bestSplitValue = 0.0;
		double maxGain = 0.0;
		//int prevTarget = dataset[0].target;
		

		// Sort the dataset based on the attribute values
		std::vector<DataPoint> sortedDataset = dataset;
		std::sort(sortedDataset.begin(), sortedDataset.end(),
				  [attributeIndex](const DataPoint& a, const DataPoint& b) {
					  return a.continuousAttributes[attributeIndex] < b.continuousAttributes[attributeIndex];
				  });
		int prevTarget = sortedDataset[0].target;
		double currentSplitValue = sortedDataset[0].continuousAttributes[attributeIndex];
		// Calculate the class counts for the entire dataset
		std::unordered_map<int, int> classCounts;
		for (const auto& data : dataset) {
			++classCounts[data.target];
		}

		// Calculate the entropy of the entire dataset
		int totalSize = dataset.size();
		double entropy = calculateEntropy2(classCounts, totalSize);

		// Iterate over the sorted dataset to find the best split value
		for (const auto& data : sortedDataset) {
			if (data.target != prevTarget) {
				// Calculate the class counts and entropy for the left and right subsets
				std::unordered_map<int, int> leftClassCounts, rightClassCounts;
				int leftSize = &data - &sortedDataset[0];
				int rightSize = totalSize - leftSize;

				for (int i = 0; i < leftSize; ++i) {
					++leftClassCounts[sortedDataset[i].target];
				}
				for (int i = leftSize; i < totalSize; ++i) {
					++rightClassCounts[sortedDataset[i].target];
				}

				double leftEntropy = calculateEntropy2(leftClassCounts, leftSize);
				double rightEntropy = calculateEntropy2(rightClassCounts, rightSize);

				// Calculate the gain of the split
				double leftWeight = static_cast<double>(leftSize) / totalSize;
				double rightWeight = static_cast<double>(rightSize) / totalSize;
				double gain = entropy - (leftWeight * leftEntropy + rightWeight * rightEntropy);

				// Update the best split value if the gain is higher
				if (gain > maxGain) {
					maxGain = gain;
					bestSplitValue = (currentSplitValue + data.continuousAttributes[attributeIndex]) / 2;
				}
				
				prevTarget = data.target;
			}

			currentSplitValue = data.continuousAttributes[attributeIndex];
		}

		return bestSplitValue;
	}


	double calculateSplitInfo(const std::vector<DataPoint>& dataset, int attributeIndex, double splitValue, bool isCategorical) {
		// Check if the attribute is categorical
		if (isCategorical) {
			// Initialize variables
			int totalSize = dataset.size();
			std::unordered_map<int, int> valueCounts;

			// Count the number of instances for each value of the categorical attribute
			for (const auto& data : dataset) {
				int value = data.categoricalAttributes[attributeIndex];
				++valueCounts[value];
			}

			// Calculate the splitInfo for categorical attributes
			double splitInfo = 0.0;
			for (const auto& entry : valueCounts) {
				double ratio = static_cast<double>(entry.second) / totalSize;
				splitInfo -= ratio * log2(ratio);
			}
			return splitInfo;
		} else {
			// Calculate the splitInfo for continuous attributes
			// Initialize variables
			int totalSize = dataset.size();
			int leftSize = 0;

			// Count the number of instances in the left subset
			for (const auto& data : dataset) {
				if (data.continuousAttributes[attributeIndex] <= splitValue) {
					++leftSize;
				}
			}

			// Calculate the ratio of the left subset size to the total size
			double leftRatio = static_cast<double>(leftSize) / totalSize;
			double rightRatio = 1.0 - leftRatio;

			
			double splitInfo = -(leftRatio * log2(leftRatio) + rightRatio * log2(rightRatio));
			return splitInfo;
		}
	}

	double calculateGainRatio(const vector<DataPoint>& dataset, int attributeIndex, const vector<bool>& usedAttributes,bool categorical) {
		double splitValue;
		double gain = calculateGain(dataset, attributeIndex,categorical, splitValue);
		double splitInfo = calculateSplitInfo(dataset, attributeIndex, splitValue, categorical);

		if (splitInfo == 0.0) {
			return 0.0;
		}

		double gainRatio = gain / splitInfo;
		return gainRatio;
	}
	
	int selectAttributeToSplit(const vector<DataPoint>& dataset, const vector<bool>& usedAttributes, bool &termination) {
		int numAttributes = usedAttributes.size();
		double averageGainRatio = 0.0;
		double averageinfogain =0.0;
		double totalGainRatio = 0.0;
		double totalinfogain = 0.0;
		int count = 0;
		bool categorical = false;
		vector<pair<double,int>> infogains;
		vector<pair<double, int>> gainRatios;  // Store gain ratios for reuse
		double splitValue;
		for (int i = 0; i < numAttributes; ++i) {
			categorical = false;
			if (i < dataset[0].categoricalAttributes.size()) {
				categorical = true;
			}
			if (!usedAttributes[i]) {
				double infogain =  calculateGain(dataset, i, categorical, splitValue);
				double gainRatio = calculateGainRatio(dataset, i, usedAttributes, categorical);
				totalGainRatio += gainRatio;
				totalinfogain += infogain;
				count++;
				
				gainRatios.emplace_back(gainRatio,i);  // Store gain ratio
				infogains.emplace_back(infogain,i);
			}
		}

		//averageGainRatio = totalGainRatio / count;
		averageinfogain = totalinfogain / count;
		vector<pair<double, int>> attributeProbabilities;
						
		for (int i = 0; i < gainRatios.size(); ++i) {
			double probability = (infogains[i].first >= averageinfogain) ? gainRatios[i].first : 0.0;
			attributeProbabilities.emplace_back(probability, gainRatios[i].second);

		}

		// Randomly select an attribute based on the probabilities

		vector<double> weights;
		int termCount = 0;//counter to check if there is no attribute with IG > 0.01
		for (const auto& pair : attributeProbabilities) {
			if(pair.first<0.01){
				termCount++;
			}
			weights.push_back(pair.first);
		}
		if(termCount==weights.size()){
			termination = true;
		}
		// Randomly select an attribute based on the probabilities
		random_device rd;
		mt19937 gen(rd());
		discrete_distribution<int> distribution(weights.begin(), weights.end());


		return attributeProbabilities[distribution(gen)].second;



	}
	vector<int> getMostCommonClass(const std::vector<DataPoint>& dataset) {
		std::unordered_map<int, int> classCounts;

		// Count the occurrences of each class label
		for (const auto& data : dataset) {
			int classLabel = data.target;  
			classCounts[classLabel]++;
		}

		// Find the most common class label
		int mostCommonClass = -1;
		int maxCount = 0;
		for (const auto& pair : classCounts) {
			if (pair.second > maxCount) {
				mostCommonClass = pair.first;
				maxCount = pair.second;
			}
		}
		vector<int> mcc = {mostCommonClass,maxCount};
		return mcc;
	}
	int missclassified(const vector<DataPoint>& dataset){
		vector<int> mcc = getMostCommonClass(dataset);
		int missclassifiedInstances = dataset.size() - mcc[1];
		return missclassifiedInstances;
	}
	Node* createDecisionTree(const vector<DataPoint>& dataset,vector<bool>& usedAttributes,const int minimumInstances) {
		Node* node = new Node();
		//cout<<dataset.size()<<endl;
		int numAttributes = usedAttributes.size();
		int numExamples = dataset.size();
		int targetClass = dataset[0].target;

		bool sameClass = true;
		bool termination = false;

		for (int i = 1; i < numExamples-1; ++i) {
			if (dataset[i].target != targetClass) {
				sameClass = false;
				break;
			}
		}

		if (sameClass) {
			node->leafNode = true;
			node->leafValue = targetClass;
			return node;
		}

		if (all_of(usedAttributes.begin(), usedAttributes.end(), [](bool used) { return used; })) {
			node->leafNode = true;
			vector<int> mcc = getMostCommonClass(dataset);
			node->leafValue = mcc[0];
			return node;
		}

		int attributeIndex = selectAttributeToSplit(dataset, usedAttributes, termination);
		node->attributeIndex = attributeIndex;
		usedAttributes[attributeIndex] = true;

		if (dataset[0].categoricalAttributes.size() > attributeIndex) {
			map<int, vector<DataPoint>> attributeValues;

			for (const auto& data : dataset) {
				attributeValues[data.categoricalAttributes[attributeIndex]].push_back(data);
			}
			//TERMINATION CRITERIA : The instances of the node’s subset are less than the minimum instances
			for (const auto& attr : attributeValues) {
				if(attr.second.size() < minimumInstances ){
					//cout<<"into minimum instances criteria minimumInstances = "<<minimumInstances<<" and child node dataset instances = "<<attr.second.size()<<endl;
					node->leafNode = true;
					vector<int> mcc = getMostCommonClass(dataset);
					node->leafValue = mcc[0];
					return node;
				}
			}

			//TERMINATION CRITERIA : The missclassified instances of the father are less or equal to children
			int childrenMissclassifiedInstances=0;
			for (const auto& entry : attributeValues) {
				childrenMissclassifiedInstances += missclassified(entry.second);
			}
			if(childrenMissclassifiedInstances >= missclassified(dataset)){
				node->leafNode = true;
				vector<int> mcc = getMostCommonClass(dataset);
				node->leafValue = mcc[0];
				return node;
			}
			//Placing the nodes as children
			for (const auto& entry : attributeValues) {
				Node* childNode = createDecisionTree(entry.second, usedAttributes,minimumInstances);
				childNode->leafValue = entry.first;
				node->children.push_back(childNode);
			}
		} else if (dataset[0].continuousAttributes.size() + dataset[0].categoricalAttributes.size() > attributeIndex) {
			double splitValue = getBestSplitValue(dataset, attributeIndex);
			node->splitValue = splitValue;

			vector<DataPoint> leftSubset;
			vector<DataPoint> rightSubset;

			for (const auto& data : dataset) {
				if (data.continuousAttributes[attributeIndex] <= splitValue) {
					leftSubset.push_back(data);
				} else {
					rightSubset.push_back(data);
				}
			}
			//TERMINATION CRITERIA : The instances of the nodes subset are less than the minimum instances
			if(leftSubset.size() < minimumInstances || rightSubset.size() < minimumInstances){
				//cout<<"into minimum instances criteria minimumInstances = "<<minimumInstances<<" and right child dataset instances = "<<rightSubset.size()<<"     and left child dataset instances = "<<leftSubset.size()<<endl;
				node->leafNode = true;
				vector<int> mcc = getMostCommonClass(dataset);
				node->leafValue = mcc[0];
				return node;
			}


			//TERMINATION CRITERIA : The missclassified instances of the father are less or equal to children

			int childrenMissclassifiedInstances=0;
			childrenMissclassifiedInstances += missclassified(leftSubset);
			childrenMissclassifiedInstances += missclassified(rightSubset);


			if(childrenMissclassifiedInstances >= missclassified(dataset)){
				node->leafNode = true;
				vector<int> mcc = getMostCommonClass(dataset);
				node->leafValue = mcc[0];
				return node;
			}


			//Placing the nodes as children
			Node* leftChildNode = createDecisionTree(leftSubset, usedAttributes,minimumInstances);
			leftChildNode->splitValue = splitValue;
			node->children.push_back(leftChildNode);

			Node* rightChildNode = createDecisionTree(rightSubset, usedAttributes,minimumInstances);
			rightChildNode->splitValue = splitValue;
			node->children.push_back(rightChildNode);
		}

		usedAttributes[attributeIndex] = false;

		return node;
	}

	int predict(const Node* node, const DataPoint& data) {
		if (node->leafNode) {
			return node->leafValue;
		}

		int attributeIndex = node->attributeIndex;
		double splitValue = node->splitValue;

		if (data.categoricalAttributes.size() > attributeIndex) {
			for (const auto& child : node->children) {
				if (child->leafValue == data.categoricalAttributes[attributeIndex]) {
					return predict(child, data);
				}
			}
		} else if (data.continuousAttributes.size() + data.categoricalAttributes.size() > attributeIndex) {
			if (data.continuousAttributes[attributeIndex] <= splitValue) {
				return predict(node->children[0], data);
			} else {
				return predict(node->children[1], data);
			}
		}

		return -1;  // Error
	}
	void printsf(const Node* node){
		if (node->leafNode) {
			cout<< "Leaf node : value = "<<node->leafValue<<endl;
		}
		else{
			cout<<"Attribute = "<<node->attributeIndex<<" and value = "<<node->leafValue<<"  if categorical || "<<node->splitValue<<"  if continuous"<<endl;
			for (const auto& child : node->children) {
					printsf(child);
			
			}
		}
	}
	int main(int argc, char** argv) {
		srand( (unsigned)time( NULL ) );
		vector<DataPoint> dataset;
		if (argc < 2)
		  return(EXIT_FAILURE);
		ifstream in(argv[1]);
		if (!in)
		  return(EXIT_FAILURE);

		loadCSV(in, dataset);
		cout<<"end of loading"<<endl;


		random_device rd;
		mt19937 g(rd());
		shuffle(dataset.begin(), dataset.end(), g);


		vector<Node*> forest;
		vector<DataPoint> testDataset,trainDataset1,trainDataset;
		for (const auto& data : dataset) {
			double chance = (double) rand()/RAND_MAX;
			if( chance > 0.8){
				testDataset.push_back(data);
			}else{
				trainDataset.push_back(data);
			}
		}
		std::unordered_map<int, int> classCounts;
		// Count the occurrences of each class label
		for (const auto& data : trainDataset) {
			int classLabel = data.target;  
			classCounts[classLabel]++;
		}
		//Calculate minimum instances
		int minimumInstances = floor(0.1 * trainDataset.size()/classCounts.size());
		cout<<"TEST DATASET SIZE "<<testDataset.size()<<" :: TRAIN DATASET SIZE "<<trainDataset.size()<<endl;

		for(int i=0; i<100; i++){
			//Splitting the dataset into train and test data 80/20
			/*
			for (const auto& data : dataset) {
				double chance = (double) rand()/RAND_MAX;
				if( chance > 0.2){
					trainDataset1.push_back(data);
				}
			}
			*/
			//vector<bool> usedAttributes(trainDataset1[0].categoricalAttributes.size() + trainDataset1[0].continuousAttributes.size(), false);
			
			vector<bool> usedAttributes(trainDataset[0].categoricalAttributes.size() + trainDataset[0].continuousAttributes.size(), false);
			//Node* root = createDecisionTree(trainDataset1, usedAttributes,minimumInstances);
			Node* root = createDecisionTree(trainDataset, usedAttributes,minimumInstances);
			
			forest.push_back(root);
			cout<<"Tree Pushed"<<endl;
			cout<<forest.size()<<endl;
			//trainDataset1.clear();
			//trainDataset.clear();
		}

		// Test the decision forest

		std::unordered_map<int, int> votes;

		int correctCount =0;
		for (const auto& testPoint : testDataset) {
			for( const auto& tree : forest){
				int vote = predict(tree,testPoint);
				votes[vote]++;
			}
			int max = -1;
			int result = -1;
			for(const auto& vote : votes){
				if(vote.second>max){
					max = vote.second;
					result = vote.first;
				}
			}
			if(result == testPoint.target){
				correctCount++;
			}
			votes.clear();
		}

		double accuracy;
		accuracy = (double)correctCount*100.00/(double)testDataset.size();
		cout<<"Accuracy = "<<accuracy<<endl;
		cout << "Correct count  = " << correctCount<< endl;

		// Clean up
		forest.clear();

		return 0;
	}

/*******************************************************************************
* Copyright (c) 2012-2014, The Microsystems Design Labratory (MDL)
* Department of Computer Science and Engineering, The Pennsylvania State University
* All rights reserved.
* 
* This source code is part of NVMain - A cycle accurate timing, bit accurate
* energy simulator for both volatile (e.g., DRAM) and non-volatile memory
* (e.g., PCRAM). The source code is free and you can redistribute and/or
* modify it by providing that the following conditions are met:
* 
*  1) Redistributions of source code must retain the above copyright notice,
*     this list of conditions and the following disclaimer.
* 
*  2) Redistributions in binary form must reproduce the above copyright notice,
*     this list of conditions and the following disclaimer in the documentation
*     and/or other materials provided with the distribution.
* 
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
* ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
* FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
* DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
* SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
* CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
* OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
* OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
* 
* Author list: 
*   Matt Poremba    ( Email: mrp5060 at psu dot edu 
*                     Website: http://www.cse.psu.edu/~poremba/ )
*******************************************************************************/

#ifndef __FRFCFS_H__
#define __FRFCFS_H__

#include "src/MemoryController.h"
#include <deque>

//EDFPCscheme
#define SAMPLECOUNT 128
#define FPCCOUNT 3
#define BDICOUNT 8
#define DYNAMICWORDSIZE 8 //chars

namespace NVM {

class FRFCFS : public MemoryController
{
  public:
    FRFCFS( );
    ~FRFCFS( );

    bool IssueCommand( NVMainRequest *req );
    bool IsIssuable( NVMainRequest *request, FailReason *fail = NULL );
    bool RequestComplete( NVMainRequest * request );

    void SetConfig( Config *conf, bool createChildren = true );

    void Cycle( ncycle_t steps );

    void RegisterStats( );
    void CalculateStats( );

  private:
    NVMTransactionQueue *memQueue;

    /* Cached Configuration Variables*/
    uint64_t queueSize;

    /* Stats */
    uint64_t measuredLatencies, measuredQueueLatencies, measuredTotalLatencies;
    double averageLatency, averageQueueLatency, averageTotalLatency;
    uint64_t mem_reads, mem_writes;
    uint64_t rb_hits;
    uint64_t rb_miss;
    uint64_t starvation_precharges;
    uint64_t cpu_insts;
    uint64_t write_pauses;
	//EDFPCscheme
    
    uint64_t bit_write;
    uint64_t bit_write_before;
    double compress_ratio;
    
    uint64_t my_llabs ( int64_t x );
    uint64_t my_abs ( int x );
    bool GeneralCompress (NVMainRequest *request, uint64_t compress);
    bool Encoder (NVMainRequest *request, bool flag);
    bool GeneralEncoder (NVMainRequest *request);
    uint64_t GetChanges (NVMainRequest *request, uint32_t MLCLevels, bool DCWFlag);
    uint64_t * convertByte2Word (NVMainRequest *request, bool flag, uint64_t size, uint64_t step);//flag: false-olddata true-newdata
    bool Word2Byte (NVMainRequest *request, bool flag, uint64_t size, uint64_t comSize, uint64_t *words, uint64_t *wordPos);//flag: false-olddata true-newdata
    
    bool FPCCompress(NVMainRequest *request, uint64_t size, bool flag );
    bool BDICompress (NVMainRequest *request, uint64_t _blockSize, bool flag );
    bool DFPCCompress(NVMainRequest *request, uint64_t _blockSize);
    
    
    bool isZeroPackable ( uint64_t * values, uint64_t size);
    bool isSameValuePackable ( uint64_t * values, uint64_t size);
    uint64_t multBaseCompression ( uint64_t * values, uint64_t size, uint64_t blimit, uint64_t bsize, uint64_t *currWords, uint64_t *currWordPos, uint64_t &pos);
    
    bool StaticCompress(NVMainRequest *request, uint64_t size, bool flag );
    uint64_t FPCIdentify (NVMainRequest *request, uint64_t size);
	uint64_t BDIIdentify (NVMainRequest *request, uint64_t _blockSize);
	uint64_t Sample (NVMainRequest *request, uint64_t _blockSize);
    uint64_t ExtractPattern();
    void HeapSort(uint64_t pattern_array[], uint64_t bytes_array[],int length, int topk);
    void HeapAdjust(uint64_t pattern_array[], uint64_t bytes_array[],int pos,int nLength);
	bool DynamicCompress(NVMainRequest *request, uint64_t size, bool flag );
    
    uint64_t DynamicFPCCompress(NVMainRequest *request, uint64_t size, bool flag );
    uint64_t DynamicBDICompress(NVMainRequest *request, uint64_t _blockSize, bool flag );
    
	uint64_t BDIpatterncounter;
    bool sample_flag;
    //uint32_t pattern_num;
	int pattern_num;
	uint64_t granularities;
    double threshold_factor;
	uint64_t FPCCounter[FPCCOUNT];
	uint64_t BDICounter[BDICOUNT];
	uint64_t SampleCounter[SAMPLECOUNT];
    uint64_t masks[FPCCOUNT+SAMPLECOUNT/DYNAMICWORDSIZE];
    int compressibleChars[FPCCOUNT+SAMPLECOUNT/DYNAMICWORDSIZE];
    bool special_pattern_flag[1+BDICOUNT];
    int mask_pos;
    
    bool encodeFlag;
    uint64_t compressIndex;
    //huffmanfpc
    bool HFPCCompress(NVMainRequest *request, uint64_t size, bool flag);



};
struct HuffmanNode {
    HuffmanNode(char k, uint64_t w) : key(k), weight(w) ,
                                 left(nullptr),
                                 right(nullptr)
    {}
    HuffmanNode(uint64_t w) : key('\0'), weight(w) ,
                                 left(nullptr),
                                 right(nullptr)
    {}
    char key;
    uint64_t weight;
    HuffmanNode * left;
    HuffmanNode * right;
};

class ComHuffmanNode {
public:
    bool operator()(HuffmanNode * e1, HuffmanNode * e2) {
        return e1->weight > e2->weight;
    }
};

class HuffmanTree {
public:
    HuffmanTree() {
        DecodeTree = nullptr;
    }
    ~HuffmanTree() {
        ClearDecodeTree();
    }

    void Input(const map<char, uint64_t> & mapCh) {
        vector<HuffmanNode*> vecHufNode;
        for (auto itr=mapCh.begin(); itr!=mapCh.end(); ++itr) {
            vecHufNode.push_back(new HuffmanNode(itr->first, itr->second));
        }

        make_heap(vecHufNode.begin(), vecHufNode.end(), ComHuffmanNode());

        while (vecHufNode.size() > 1) {
            HuffmanNode * right = vecHufNode.front();
            pop_heap(vecHufNode.begin(), vecHufNode.end(), ComHuffmanNode());
            vecHufNode.pop_back();

            HuffmanNode * left = vecHufNode.front();
            pop_heap(vecHufNode.begin(), vecHufNode.end(), ComHuffmanNode());
            vecHufNode.pop_back();

            HuffmanNode * parent = new HuffmanNode(left->weight + right->weight);
            parent->left = left;
            parent->right = right;

            vecHufNode.push_back(parent);
            push_heap(vecHufNode.begin(), vecHufNode.end(), ComHuffmanNode());
        }

        if (!vecHufNode.empty()) {
            DecodeTree = vecHufNode.front();
        }

        veccode.resize(std::numeric_limits<char>().max());

        string code;

        BuildCode(DecodeTree, code);

    }

    void ClearDecodeTree(HuffmanNode * pNode) {
        if (pNode == nullptr) return;

        ClearDecodeTree(pNode->left);
        ClearDecodeTree(pNode->right);
        delete pNode;
    }

    void ClearDecodeTree() {
        ClearDecodeTree(DecodeTree);
        DecodeTree = nullptr;
    }

    void BuildCode(HuffmanNode * pNode, string & code) {
        if (pNode->left == NULL) {
            veccode[pNode->key] = code;
            return ;
        }

        code.push_back('0');
        BuildCode(pNode->left, code);
        code.pop_back();
        code.push_back('1');
        BuildCode(pNode->right, code);
        code.pop_back();
    }

    string Decode(const string & strB) {
        string strC;

        HuffmanNode * pNode = DecodeTree;
        for (unsigned int i=0; i<strB.size(); ++i) {
            if (strB[i] == '0') {
                pNode = pNode->left;
            } else {
                pNode = pNode->right;
            }

            if (pNode->left == NULL) {
                strC.push_back(pNode->key);
                pNode = DecodeTree;
            }
        }

        return strC;
    }

    string GetCode(const string & strA) {
        string strB;
        for (unsigned int i=0; i<strA.size(); ++i) {
            strB += veccode[strA[i]];
        }

        return strB;
    }

    vector<string> veccode;
private:
    HuffmanNode * DecodeTree;
    
};


};

#endif

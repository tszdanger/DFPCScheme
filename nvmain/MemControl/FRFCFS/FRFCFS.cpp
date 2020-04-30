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

#include "MemControl/FRFCFS/FRFCFS.h"
#include "src/EventQueue.h"
#include "include/NVMainRequest.h"
#ifndef TRACE
#ifdef GEM5
  #include "SimInterface/Gem5Interface/Gem5Interface.h"
  #include "base/statistics.hh"
  #include "base/types.hh"
  #include "sim/core.hh"
  #include "sim/stat_control.hh"
#endif
#endif
#include <iostream>
#include <set>
#include <assert.h>
#include <stdio.h>
#include <vector>
#include <string>
#include <map>
#include <algorithm>
using namespace std;
// struct HuffmanNode {
//     HuffmanNode(char k, uint64_t w) : key(k), weight(w) ,
//                                  left(nullptr),
//                                  right(nullptr)
//     {}
//     HuffmanNode(uint64_t w) : key('\0'), weight(w) ,
//                                  left(nullptr),
//                                  right(nullptr)
//     {}
//     char key;
//     uint64_t weight;
//     HuffmanNode * left;
//     HuffmanNode * right;
// };

// class ComHuffmanNode {
// public:
//     bool operator()(HuffmanNode * e1, HuffmanNode * e2) {
//         return e1->weight > e2->weight;
//     }
// };


bool ComHuffmanNode::operator()(HuffmanNode * e1, HuffmanNode * e2) {
        return e1->weight > e2->weight;
    }

HuffmanTree::HuffmanTree(){DecodeTree = nullptr;}

HuffmanTree::~HuffmanTree() {
        ClearDecodeTree();
    }

void HuffmanTree::Input(const std::map<char, uint64_t> & mapCh) {
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

void HuffmanTree::ClearDecodeTree(HuffmanNode * pNode) {
        if (pNode == nullptr) return;

        ClearDecodeTree(pNode->left);
        ClearDecodeTree(pNode->right);
        delete pNode;
    }

void HuffmanTree::ClearDecodeTree() {
        ClearDecodeTree(DecodeTree);
        DecodeTree = nullptr;
    }

void HuffmanTree::BuildCode(HuffmanNode * pNode, std::string & code) {
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

std::string HuffmanTree::Decode(const string & strB) {
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

std::string HuffmanTree::GetCode(const string & strA) {
        string strB;
        for (unsigned int i=0; i<strA.size(); ++i) {
            strB += veccode[strA[i]];
        }

        return strB;
    }


// class HuffmanTree {
// public:
//     HuffmanTree() {
//         DecodeTree = nullptr;
//     }
//     ~HuffmanTree() {
//         ClearDecodeTree();
//     }

//     void Input(const std::map<char, uint64_t> & mapCh) {
//         vector<HuffmanNode*> vecHufNode;
//         for (auto itr=mapCh.begin(); itr!=mapCh.end(); ++itr) {
//             vecHufNode.push_back(new HuffmanNode(itr->first, itr->second));
//         }

//         make_heap(vecHufNode.begin(), vecHufNode.end(), ComHuffmanNode());

//         while (vecHufNode.size() > 1) {
//             HuffmanNode * right = vecHufNode.front();
//             pop_heap(vecHufNode.begin(), vecHufNode.end(), ComHuffmanNode());
//             vecHufNode.pop_back();

//             HuffmanNode * left = vecHufNode.front();
//             pop_heap(vecHufNode.begin(), vecHufNode.end(), ComHuffmanNode());
//             vecHufNode.pop_back();

//             HuffmanNode * parent = new HuffmanNode(left->weight + right->weight);
//             parent->left = left;
//             parent->right = right;

//             vecHufNode.push_back(parent);
//             push_heap(vecHufNode.begin(), vecHufNode.end(), ComHuffmanNode());
//         }

//         if (!vecHufNode.empty()) {
//             DecodeTree = vecHufNode.front();
//         }

//         veccode.resize(std::numeric_limits<char>().max());

//         string code;

//         BuildCode(DecodeTree, code);

//     }

//     void ClearDecodeTree(HuffmanNode * pNode) {
//         if (pNode == nullptr) return;

//         ClearDecodeTree(pNode->left);
//         ClearDecodeTree(pNode->right);
//         delete pNode;
//     }

//     void ClearDecodeTree() {
//         ClearDecodeTree(DecodeTree);
//         DecodeTree = nullptr;
//     }

//     void BuildCode(HuffmanNode * pNode, std::string & code) {
//         if (pNode->left == NULL) {
//             veccode[pNode->key] = code;
//             return ;
//         }

//         code.push_back('0');
//         BuildCode(pNode->left, code);
//         code.pop_back();
//         code.push_back('1');
//         BuildCode(pNode->right, code);
//         code.pop_back();
//     }

//     std::string Decode(const string & strB) {
//         string strC;

//         HuffmanNode * pNode = DecodeTree;
//         for (unsigned int i=0; i<strB.size(); ++i) {
//             if (strB[i] == '0') {
//                 pNode = pNode->left;
//             } else {
//                 pNode = pNode->right;
//             }

//             if (pNode->left == NULL) {
//                 strC.push_back(pNode->key);
//                 pNode = DecodeTree;
//             }
//         }

//         return strC;
//     }

//     std::string GetCode(const string & strA) {
//         string strB;
//         for (unsigned int i=0; i<strA.size(); ++i) {
//             strB += veccode[strA[i]];
//         }

//         return strB;
//     }

//     std::vector<string> veccode;
// private:
//     HuffmanNode * DecodeTree;
    
// };

using namespace NVM;

FRFCFS::FRFCFS( )
{
    std::cout << "Created a First Ready First Come First Serve memory controller!"
        << std::endl;
	
	//EDFPCscheme
    encodeFlag = false;
    compressIndex = 5;//0:DCW 1:FPC 2:BDI 3:DFPC 4:HFPC 5:DHFPC
    bit_write_before = 0;
    bit_write = 0;
    
    BDIpatterncounter = 0;
    sample_flag = true;
    pattern_num = SAMPLECOUNT/DYNAMICWORDSIZE/2;
    threshold_factor = 0.4;
	granularities = 50000;
    compress_ratio = 0.0f;
    average_compress_ratio = 0.0f;
    total_compress_time = 1;
    mask_pos = 0;
    
	for(int i=0;i<FPCCOUNT;i++)
		FPCCounter[i] = 0;
	for(int i=0;i<BDICOUNT;i++)
		BDICounter[i] = 0;
	for(int i=0;i<SAMPLECOUNT;i++)
		SampleCounter[i] = 0;
    for(int i = 0; i <= BDICOUNT; i++)
        special_pattern_flag[i] = false;
	
    queueSize = 32;
    starvationThreshold = 4;

    averageLatency = 0.0f;
    averageQueueLatency = 0.0f;
    averageTotalLatency = 0.0f;

    measuredLatencies = 0;
    measuredQueueLatencies = 0;
    measuredTotalLatencies = 0;

    mem_reads = 0;
    mem_writes = 0;

    rb_hits = 0;
    rb_miss = 0;

    write_pauses = 0;

    starvation_precharges = 0;
    
    
    //默认没有建树
    buildtree = false;

    

    psInterval = 0;

    InitQueues( 1 );

    memQueue = &(transactionQueues[0]);
}

FRFCFS::~FRFCFS( )
{
    std::cout << "FRFCFS memory controller destroyed. " << memQueue->size( ) 
              << " commands still in memory queue." << std::endl;
}

void FRFCFS::SetConfig( Config *conf, bool createChildren )
{
    if( conf->KeyExists( "StarvationThreshold" ) )
    {
        starvationThreshold = static_cast<unsigned int>( conf->GetValue( "StarvationThreshold" ) );
    }

    if( conf->KeyExists( "QueueSize" ) )
    {
        queueSize = static_cast<unsigned int>( conf->GetValue( "QueueSize" ) );
    }

    MemoryController::SetConfig( conf, createChildren );

    SetDebugName( "FRFCFS", conf );
}

void FRFCFS::RegisterStats( )
{
	
	//EDFPCscheme
    AddStat(bit_write);
    AddStat(bit_write_before);
    AddStat(compress_ratio);
    AddStat(average_compress_ratio);
	
    AddStat(mem_reads);
    AddStat(mem_writes);
    AddStat(rb_hits);
    AddStat(rb_miss);
    AddStat(starvation_precharges);
    AddStat(averageLatency);
    AddStat(averageQueueLatency);
    AddStat(averageTotalLatency);
    AddStat(measuredLatencies);
    AddStat(measuredQueueLatencies);
    AddStat(measuredTotalLatencies);
    AddStat(write_pauses);

    MemoryController::RegisterStats( );
}

bool FRFCFS::IsIssuable( NVMainRequest * /*request*/, FailReason * /*fail*/ )
{
    bool rv = true;

    /*
     *  Limit the number of commands in the queue. This will stall the caches/CPU.
     */ 
    if( memQueue->size( ) >= queueSize )
    {
        rv = false;
    }

    return rv;
}

/*
 *  This method is called whenever a new transaction from the processor issued to
 *  this memory controller / channel. All scheduling decisions should be made here.
 */
bool FRFCFS::IssueCommand( NVMainRequest *req )
{
    if( !IsIssuable( req ) )
    {
        return false;
    }

    req->arrivalCycle = GetEventQueue()->GetCurrentCycle();

    /* 
     *  Just push back the read/write. It's easier to inject dram commands than break it up
     *  here and attempt to remove them later.
     */
     //EDFPC
	if( req->type == WRITE)
	{
        //bool isCom = false;
        //bool isEncoded = false;
        
        uint64_t bitsChange = 0;
        uint64_t bitsCom = 0;
        // uint64_t i;
        uint64_t size;
        uint64_t comsize;
        // char str[20];
        comsize = size = req->data.GetSize();
        // printf("data: ");
        // for(uint64_t i = 0; i < size; i++)
        // {
        //     printf("%02X", req->data.GetByte(i));
        //     // sprintf(str, 20, "%X", req->data.GetByte(i));
        //     // puts(str);

        // }
        // printf("\n");
        
		//isCom = GeneralCompress(req, 3);
        
        
        //统计一下这些模式有多少个,仅FPC那种不包含bdi
        
        uint64_t * values1 = convertByte2Word(req, false, 64, 8);
        if( isZeroPackable( values1, 64 / 8))
        {
            pattern_ana[8]+=1;
        }
        free(values1);
        values1 = convertByte2Word(req, false, size*4, 4);
        for (int i = 0; i < size; i++)
    {
        if(values1[i]==0){
            pattern_ana[0]++;
            continue;
        }
        if(my_abs((int)(values1[i])) <= 0xFF){
            pattern_ana[1]++;
            continue;
        }
        if(my_abs((int)(values1[i])) <= 0xFFFF){
            pattern_ana[2]++;
            continue;
        }
        if(((values1[i]) & 0xFFFF) == 0 ){
            pattern_ana[3]++;
            continue;
        }
        if( my_abs((int)((values1[i]) & 0xFFFF)) <= 0xFF
             && my_abs((int)((values1[i] >> 16) & 0xFFFF)) <= 0xFF){
            pattern_ana[4]++;
            continue;
        }
        if(my_abs((int)(values1[i])) <= 0xF){
            pattern_ana[7]++;
            continue;
        }
        uint64_t byte0 = (values1[i]) & 0xFF;
        uint64_t byte1 = (values1[i] >> 8) & 0xFF;
        uint64_t byte2 = (values1[i] >> 16) & 0xFF;
        uint64_t byte3 = (values1[i] >> 24) & 0xFF;
        if(byte0 == byte1 && byte0 == byte2 && byte0 == byte3){
            pattern_ana[5]++;
            continue;
        }
        pattern_ana[6]++;
    }





        GeneralCompress(req, compressIndex);
        if(req->data.IsCompressed())
        {
            comsize = req->data.GetComSize();
        }
        compress_ratio += (size * 1.0 / comsize);
        average_compress_ratio = ((size * 1.0 / comsize)+average_compress_ratio*(total_compress_time-1))/(total_compress_time);
        total_compress_time++;
        //100000输出一下模式统计
        if(total_compress_time%100000==0)
        {
         std::cout<< "size is "<< size<<"comsize is"<< comsize << "so compress_ratio is"<< compress_ratio <<std::endl;
         std::cout<<"average_compress_ratio is"<<average_compress_ratio<<"total_compression time is"<<total_compress_time<<std::endl;
         std::cout<<"averagelatency is"<<averageLatency<<std::endl;
         
        //  std::cout<< "模式数目分别为: "<<std::endl;
        //  for(int j=0;j<9;j++){
        //     std::cout<< pattern_ana[j] <<"\t";
        //  }
         printf("\n");
        }



        /*
        if(req->data.IsCompressed())
        {
            size = req->data.GetComSize();
            printf("comdata: ");
            for(i = 0; i<size; i++)
            {
                printf("%d ", req->data.GetComByte(i));
            }
            printf("\n");
        }
        
        if(req->oldData.IsCompressed())
        {
            printf("comolddata: ");
            size = req->oldData.GetComSize();
            for(i = 0; i<size; i++)
            {
                printf("%d ", req->oldData.GetComByte(i));
            }
            printf("\n");
        }
    */
        if(encodeFlag)
        {
            //isEncoded = GeneralEncoder(req);
            GeneralEncoder(req);
            /*
            if(req->data.IsCompressed())
            {
                printf("encodedata: ");
                size = req->data.GetComSize();
                for(i = 0; i<size; i++)
                {
                    printf("%d ", req->data.GetComByte(i));
                }
                printf("\n");
            }
            if(req->oldData.IsCompressed())
            {
                printf("encodeolddata: ");
                size = req->oldData.GetComSize();
                for(i = 0; i<size; i++)
                {
                    printf("%d ", req->oldData.GetComByte(i));
                }
                printf("\n");
            }
            */
        }
        //if(isEncoded)
        //{
        //    
        //}
        
        
        bitsCom = GetChanges(req, 3, false);
        bitsChange = GetChanges(req, 3, true);
        //std::cout<<bitsChange<<std::endl;
        if(bitsChange > 64 * 8)
        {
            //error
            bitsChange = 64 * 8;
        }
        //std::cout<<"bit_write:"<<bit_write<<std::endl;
        bit_write_before += bitsCom;
		bit_write += bitsChange;
        //std::cout<<"bit_write_com:"<<bit_write<<std::endl;
	}
	 
    Enqueue( 0, req );

    if( req->type == READ )
        mem_reads++;
    else
        mem_writes++;

    /*
     *  Return whether the request could be queued. Return false if the queue is full.
     */
    return true;
}

bool FRFCFS::RequestComplete( NVMainRequest * request )
{
    if( request->type == WRITE || request->type == WRITE_PRECHARGE )
    {
        /* 
         *  Put cancelled requests at the head of the write queue
         *  like nothing ever happened.
         */
        if( request->flags & NVMainRequest::FLAG_CANCELLED 
            || request->flags & NVMainRequest::FLAG_PAUSED )
        {
            Prequeue( 0, request );

            return true;
        }
    }

    /* Only reads and writes are sent back to NVMain and checked for in the transaction queue. */
    if( request->type == READ 
        || request->type == READ_PRECHARGE 
        || request->type == WRITE 
        || request->type == WRITE_PRECHARGE )
    {
        request->status = MEM_REQUEST_COMPLETE;
        request->completionCycle = GetEventQueue()->GetCurrentCycle();

        /* Update the average latencies based on this request for READ/WRITE only. */
        averageLatency = ((averageLatency * static_cast<double>(measuredLatencies))
                           + static_cast<double>(request->completionCycle)
                           - static_cast<double>(request->issueCycle))
                       / static_cast<double>(measuredLatencies+1);
        measuredLatencies += 1;

        averageQueueLatency = ((averageQueueLatency * static_cast<double>(measuredQueueLatencies))
                                + static_cast<double>(request->issueCycle)
                                - static_cast<double>(request->arrivalCycle))
                            / static_cast<double>(measuredQueueLatencies+1);
        measuredQueueLatencies += 1;

        averageTotalLatency = ((averageTotalLatency * static_cast<double>(measuredTotalLatencies))
                                + static_cast<double>(request->completionCycle)
                                - static_cast<double>(request->arrivalCycle))
                            / static_cast<double>(measuredTotalLatencies+1);
        measuredTotalLatencies += 1;
    }

    return MemoryController::RequestComplete( request );
}

void FRFCFS::Cycle( ncycle_t steps )
{
    NVMainRequest *nextRequest = NULL;

    /* Check for starved requests BEFORE row buffer hits. */
    if( FindStarvedRequest( *memQueue, &nextRequest ) )
    {
        rb_miss++;
        starvation_precharges++;
    }
    /* Check for row buffer hits. */
    else if( FindRowBufferHit( *memQueue, &nextRequest) )
    {
        rb_hits++;
    }
    /* Check if the address is accessible through any other means. */
    else if( FindCachedAddress( *memQueue, &nextRequest ) )
    {
    }
    else if( FindWriteStalledRead( *memQueue, &nextRequest ) )
    {
        if( nextRequest != NULL )
            write_pauses++;
    }
    /* Find the oldest request that can be issued. */
    else if( FindOldestReadyRequest( *memQueue, &nextRequest ) )
    {
        rb_miss++;
    }
    /* Find requests to a bank that is closed. */
    else if( FindClosedBankRequest( *memQueue, &nextRequest ) )
    {
        rb_miss++;
    }
    else
    {
        nextRequest = NULL;
    }

    /* Issue the commands for this transaction. */
    if( nextRequest != NULL )
    {
        IssueMemoryCommands( nextRequest );
    }

    /* Issue any commands in the command queues. */
    CycleCommandQueues( );

    MemoryController::Cycle( steps );
}

void FRFCFS::CalculateStats( )
{
    MemoryController::CalculateStats( );
}



//EDFPCscheme

uint64_t FRFCFS::my_llabs ( int64_t x )
{
    uint64_t t = x >> 63;
    return (x ^ t) - t;
}

uint64_t FRFCFS::my_abs ( int x )
{
    uint64_t t = x >> 31;
    return (x ^ t) - t;
}

uint64_t * FRFCFS::convertByte2Word (NVMainRequest *request, bool flag, uint64_t size, uint64_t step)//flag: false-olddata true-newdata
{
    uint64_t * values = (uint64_t *) malloc(sizeof(uint64_t) * size/step);
    uint64_t i,j; 
    for (i = 0; i < size / step; i++) {
        values[i] = 0;    // Initialize all elements to zero.
    }
    
    for (i = 0; i < size; i += step ){
        for (j = 0; j < step; j++){
            if(flag)
                values[i / step] += (uint64_t)(request->data.GetByte(i + j) << (8*j));
            else
                values[i / step] += (uint64_t)(request->oldData.GetByte(i + j) << (8*j));
                
        }
    }
    return values;
}

uint64_t FRFCFS::GetChanges (NVMainRequest *request, uint32_t MLCLevels, bool DCWFlag)
{
    uint32_t *rawData = NULL;
    uint32_t *oldData = NULL;
    uint64_t memoryWordSize = 64*8;
    uint64_t size = 0;
    uint64_t bitsChange = 0;
    uint64_t cellsChange = 0;
    if(request->data.IsCompressed())
    {
        rawData = reinterpret_cast<uint32_t*>(request->data.comData);
        memoryWordSize = request->data.GetComSize()*8;
    }
    else
    {
        rawData = reinterpret_cast<uint32_t*>(request->data.rawData);
    }
    
    if(request->oldData.IsCompressed())
    {
        oldData = reinterpret_cast<uint32_t*>(request->oldData.comData);
    }
    else
    {
        oldData = reinterpret_cast<uint32_t*>(request->oldData.rawData);
    }
    if(!DCWFlag)
        return memoryWordSize;
    size = memoryWordSize/32;
    
    if(size > 16)
        size = 16;
    uint64_t i = 0;
    ncounter_t i_pos = 0;
    ncounter_t nums = 32 / MLCLevels;
    if(32 % MLCLevels != 0)
        nums++;
    uint32_t word;
    uint32_t oldWord;
    uint32_t mask = 0x00000001;
    uint32_t cellData;
    uint32_t oldCellData;
    if(MLCLevels == 1)
        mask = 0x00000001;
    else if(MLCLevels == 2)
        mask = 0x00000003;
    else if(MLCLevels == 3)
        mask = 0x00000007;
    for(i = 0; i < size; i++)
    {
        word = rawData[i];
        oldWord = oldData[i];
        for(i_pos = 0; i_pos < nums; i_pos++)
        {
            cellData = word & mask;
            oldCellData = oldWord & mask;
            if(cellData != oldCellData)
            {
                if(i_pos == nums - 1)
                    bitsChange += 2;
                else
                    cellsChange++;
            }
            word = word >> MLCLevels; 
            oldWord = oldWord >> MLCLevels; 
        }
    }
    size = memoryWordSize % 32;
    if(size != 0)
    {
        word = rawData[i];
        oldWord = oldData[i];
        ncounter_t nums_r = size / MLCLevels;
        if(size % MLCLevels != 0)
            nums_r++;
        for(i_pos = 0; i_pos < (nums-nums_r); i_pos++)
        {
            word = word >> MLCLevels; 
            oldWord = oldWord >> MLCLevels; 
        }
        for(; i_pos < nums; i_pos++)
        {
            cellData = word & mask;
            oldCellData = oldWord & mask;
            if(cellData != oldCellData)
            {
                if(i_pos == nums - 1)
                    bitsChange += 2;
                else
                    cellsChange++;
            }
            word = word >> MLCLevels; 
            oldWord = oldWord >> MLCLevels; 
        }
    }
    bitsChange += cellsChange * MLCLevels;
    return bitsChange;
}

bool FRFCFS::Encoder (NVMainRequest *request, bool flag)
{
    bool resFlag = false;
    uint64_t size;
    bool rawDataArray[64*8];
    bool encodedDataArray[64*8];
    uint64_t i, j, iMax, jMax;
    bool isHalf;
    uint8_t dataByte;
    
    uint64_t comSize = 0;
    
    for(i = 0; i< 64*8; i++)
    {
        rawDataArray[i] = false;
        encodedDataArray[i] = false;
    }
    
    if(flag)
    {
        size = request->data.GetComSize();
        isHalf = request->data.IsHalf();
    }
    else
    {
        size = request->oldData.GetComSize();
        isHalf = request->oldData.IsHalf();
    }
        
    /*
    cells = size * 4 / 3;
    if(size * 4 % 3 != 0)
        cells++;
    if(cells <= 64)
        //encode
    */
    if(size <= 48)
    {
        resFlag = true;
        iMax = size;
        jMax = 8;
        for(i = 0; i < iMax; i++)
        {
            if(flag)
                dataByte = request->data.GetComByte(i);
            else
                dataByte = request->oldData.GetComByte(i);
            if(i == iMax - 1 && isHalf)
                jMax = 4;
            for(j = 0; j<jMax; j++)
            {
                if(((dataByte >> (7-j)) & 0x1) != 0)
                    rawDataArray[i*8+j] = true;
                else
                    rawDataArray[i*8+j] = false;
            }
        }
        iMax = size * 4;
        if(isHalf)
            iMax = iMax - 2;
        //*8/2 or (*8-4)/2
        
        for(i = 0, j = 0; i < iMax; i=i+2)
        {
            if(!rawDataArray[i*2] && !rawDataArray[i*2+1])//00
            {
                //encodedCell = 0
                encodedDataArray[j++] = 0;
                encodedDataArray[j++] = 0;
                encodedDataArray[j++] = 0;
            }
            else if(!rawDataArray[i*2] && rawDataArray[i*2+1])//01
            {
                //encodedCell = 1
                encodedDataArray[j++] = 0;
                encodedDataArray[j++] = 0;
                encodedDataArray[j++] = 1;
            }
            else if(rawDataArray[i*2] && !rawDataArray[i*2+1])//10
            {
                //encodedCell = 2
                encodedDataArray[j++] = 1;
                encodedDataArray[j++] = 1;
                encodedDataArray[j++] = 0;
            }
            else//11
            {
                //encodedCell = 3
                encodedDataArray[j++] = 1;
                encodedDataArray[j++] = 1;
                encodedDataArray[j++] = 1;
            } 
        }
        comSize = j / 8;
        if(j%8 != 0)
            comSize++;
        if(flag)
            request->data.SetComSize(comSize);
        iMax = comSize;
        jMax = 8;
        for(i = 0; i < iMax; i++)
        {
            dataByte = 0;
            for(j = 0; j<jMax; j++)
            {
                if(encodedDataArray[i*jMax+j])
                {
                    dataByte += (1 << (7-i));
                }
            }
            if(flag)
                request->data.SetComByte(i, dataByte);
            else
                request->oldData.SetComByte(i, dataByte);
        }
    }
    return resFlag;
}

bool FRFCFS::GeneralEncoder (NVMainRequest *request)
{
    if(request->data.IsCompressed())
    {
        Encoder(request, true);
    }
    if(request->oldData.IsCompressed())
    {
        Encoder(request, false);
    }
    return true;
}

bool FRFCFS::GeneralCompress (NVMainRequest *request, uint64_t compress)
{   // compress is the actual compression algorithm
    bool resFlag = false;
    uint64_t _blockSize = request->data.GetSize();//64
    switch (compress)
    {
        case 0:
            //
            break;
        case 1:	
            //FPC
            resFlag = FPCCompress(request, _blockSize/4, true);
            FPCCompress(request, _blockSize/4, false);
            break;
        case 2:
            //BDI
			resFlag = BDICompress(request, _blockSize, true);
            BDICompress(request, _blockSize, false);
            break;
        case 3:
            //DFPC
			resFlag = DFPCCompress(request, _blockSize);
            break;
        case 4:
            resFlag = HFPCCompress(request,_blockSize/4,true);
            HFPCCompress(request,_blockSize/4,false);
            break;
        case 5:
            resFlag = HDFPCCompress(request, _blockSize);
            break;
        default:
            break;
    }
    return resFlag;
}

bool FRFCFS::Word2Byte (NVMainRequest *request, bool flag, uint64_t size, uint64_t comSize, uint64_t *words, uint64_t *wordPos)//flag: false-olddata true-newdata
{
    uint64_t i,j;
    uint8_t dataChar = 0;
    uint8_t dataByte = 0;
    bool dataFlag = false;//false--low true--high
    uint64_t bytePos = 0;
    
    if(flag)//compressible newdata
    {
        request->data.SetComSize(comSize);
    }
    else
    {
        request->oldData.SetComSize(64);
    }
    
    for (i = 0; i < size; i++)
    {
        for(j = wordPos[i]; j>0; j--)
        {
            dataChar = (words[i] >> ((j-1)*4)) & 0xF;
            if(dataFlag)
            {
                dataByte = dataByte | dataChar;
                if(flag)
                {
                    request->data.SetComByte(bytePos, dataByte);
                }
                else
                {
                    request->oldData.SetComByte(bytePos, dataByte);
                }
                bytePos++;
                //std::cout<<dataByte<<" ";
            }
            else
            {
                dataByte = dataChar << 4;
            }
            dataFlag = !dataFlag;
        }
    }
    if(dataFlag)
    {
        if(flag)
        {
            request->data.SetComByte(bytePos, dataByte);
            request->data.SetHalfFlag(true);
        }
        else
        {
            request->oldData.SetComByte(bytePos, dataByte);
            request->oldData.SetHalfFlag(true);
        }
        bytePos++;
        //std::cout<<dataByte;
    }
    //std::cout<<std::endl;
    
    if(flag && (bytePos<=comSize))//compressible newdata
    {
        // std::cout << "yes! bytePos is"<<bytePos<<"and comsize is "<<comSize <<std::endl;
        request->data.SetComSize(bytePos);
    }
    //dataFlag = !dataFlag;
    
    return true;
}
//huffman fpc static version1
bool FRFCFS::Word2ByteHuff (NVMainRequest *request, bool flag, uint64_t size, uint64_t comSize, uint64_t *words, uint64_t *wordPos,int huffbit)//flag: false-olddata true-newdata
{
    uint64_t i,j;
    uint8_t dataChar = 0;
    uint8_t dataByte = 0;
    bool dataFlag = false;//false--low true--high
    uint64_t bytePos = 0;
    
    if(flag)//compressible newdata
    {
        request->data.SetComSize(comSize);
    }
    else
    {
        request->oldData.SetComSize(64);
    }
    
    for (i = 0; i < size; i++)
    {
        for(j = wordPos[i]; j>0; j--)
        {
            dataChar = (words[i] >> ((j-1)*4)) & 0xF;
            if(dataFlag)
            {
                dataByte = dataByte | dataChar;
                if(flag)
                {
                    request->data.SetComByte(bytePos, dataByte);
                }
                else
                {
                    request->oldData.SetComByte(bytePos, dataByte);
                }
                bytePos++;
                //std::cout<<dataByte<<" ";
            }
            else
            {
                dataByte = dataChar << 4;
            }
            dataFlag = !dataFlag;
        }
    }
    if(dataFlag)
    {
        if(flag)
        {
            request->data.SetComByte(bytePos, dataByte);
            request->data.SetHalfFlag(true);
        }
        else
        {
            request->oldData.SetComByte(bytePos, dataByte);
            request->oldData.SetHalfFlag(true);
        }
        bytePos++;
        //std::cout<<dataByte;
    }
    //std::cout<<std::endl;
    
    if(flag && ((bytePos+huffbit/4)<=comSize))//compressible newdata
    {
        // std::cout << "yes! bytePos+huffbit/4"<<bytePos+huffbit/4<<"and comsize is "<<comSize <<std::endl;
        request->data.SetComSize((bytePos+huffbit/4)+1);
    }
    //dataFlag = !dataFlag;
    
    return true;
}

bool FRFCFS::PUREHFPCCompress(NVMainRequest *request, uint64_t size, bool flag){
    //GetHuffCode,对64B的cacheline ，16*32bits 每32bits 采样一下，总共采样16个看pattern频率
    // 输出 all zero/8bits符号扩展/16bits符号扩展/32bits/00ab00cd类/重复类/不可压缩类 
    //here outside we still trans size/4 ,so here size*4 means size outside
    //then values[i] is a 32bits thing
    // printf("进入了PUREHFPCCompress\n");

    uint64_t * values = convertByte2Word(request, flag, size*4, 8);
    uint64_t i;
    
    uint64_t words[16];
    uint64_t wordPos[16]; //0~8 chars
    uint64_t comSize = 0;
    bool comFlag = false;
    
    if( isZeroPackable( values, size*4 / 8))
    {
        // 000
        // 000 静态全0压缩
        std::cout<<"64 all zeros" <<std::endl;
        words[0] = 0;
        wordPos[0] = 1;
        comFlag = true;
        comSize = 1;
        free(values);
        values = NULL;
        Word2Byte(request, flag, comSize, comSize, words, wordPos);
        return comFlag;
    }
    free(values);
    values = convertByte2Word(request, flag, size*4, 4);
    int huffbit=0;
    char name[8] = {'a','b','c','d','e','f','g','h'};

    for (i=0;i<size;i++){
        //这里word[i]先只存实际值,用huffbit来存使用的bits数
        //00 -> all zero
        if(values[i] == 0){
            // words[i] = values[i] + 0x0;
            words[i] = values[i] + strtol(huffTree1.veccode[name[0]].c_str(),NULL,2);
            // wordPos[i] = 1;
            wordPos[i] = 0;
            comSize += wordPos[i];
            // huffbit +=2;
            huffbit += huffTree1.veccode[name[0]].size();
            continue;
        }

        // 01 8bits符号扩展
        if(my_abs((int)(values[i])) <= 0xFF){
            words[i] = my_abs((int)(values[i])) + (strtol(huffTree1.veccode[name[1]].c_str(),NULL,2)<<8);
            wordPos[i] = 3;
            comSize += wordPos[i];
            // huffbit +=2;
            huffbit += huffTree1.veccode[name[1]].size();
            continue;
        }
        if(my_abs((int)(values[i])) <= 0xF){
            words[i] = my_abs((int)(values[i])) + (strtol(huffTree1.veccode[name[7]].c_str(),NULL,2)<<4);
            wordPos[i] = 2;
            comSize += wordPos[i];
            // huffbit +=2;
            huffbit += huffTree1.veccode[name[7]].size();
            continue;
        }
        // 110 16bits符号扩展
        if(my_abs((int)(values[i])) <= 0xFFFF){
            // words[i] = my_abs((int)(values[i])) + 0x40000;
            words[i] = my_abs((int)(values[i])) + (strtol(huffTree1.veccode[name[2]].c_str(),NULL,2)<<16);
            wordPos[i] = 4;
            comSize += wordPos[i];
            // huffbit +=3;
            huffbit += huffTree1.veccode[name[2]].size();
            continue;
        }
        //100  32bits
        if(((values[i]) & 0xFFFF) == 0 ){
            // words[i] = (values[i] >> 16) + 0x40000;
            words[i] = (values[i] >> 16) + (strtol(huffTree1.veccode[name[3]].c_str(),NULL,2)<<16);
            wordPos[i] = 4;
            comSize += wordPos[i];
            // huffbit +=3;
            huffbit += huffTree1.veccode[name[3]].size();
            continue;
        }
        //101 00ab00cd类
        if( my_abs((int)((values[i]) & 0xFFFF)) <= 0xFF
             && my_abs((int)((values[i] >> 16) & 0xFFFF)) <= 0xFF){
            words[i] = my_abs((int)((values[i] >> 8))) + my_abs((int)((values[i]) & 0xFFFF)) + (strtol(huffTree1.veccode[name[4]].c_str(),NULL,2)<<32);
            wordPos[i] = 4;
            comSize += wordPos[i];
            // huffbit +=3;
            huffbit += huffTree1.veccode[name[4]].size();
            continue;
        }
        //1110 重复类
        uint64_t byte0 = (values[i]) & 0xFF;
        uint64_t byte1 = (values[i] >> 8) & 0xFF;
        uint64_t byte2 = (values[i] >> 16) & 0xFF;
        uint64_t byte3 = (values[i] >> 24) & 0xFF;
        if(byte0 == byte1 && byte0 == byte2 && byte0 == byte3){
            words[i] = byte0 +  (strtol(huffTree1.veccode[name[5]].c_str(),NULL,2)<<8);
            wordPos[i] = 2;
            comSize += wordPos[i];
            // huffbit +=4;
            huffbit += huffTree1.veccode[name[5]].size();
            continue;
        }
        //1111
        words[i] = values[i];
        wordPos[i] = 8;
        //这里huffbit不用加因为有一bit是是否压缩的
        comSize += wordPos[i];
    }
    //这里我们把prefix加上去
    std::cout << "now huffbit is"<<huffbit <<std::endl;
    comSize += huffbit/4;
    if (comSize==0)
        comSize=1;
    if(comSize % 2 == 1)
        comSize++;
    comSize /= 2;
    if(comSize < (size*4))
    {
        comFlag = true;
    }
    
    //6 bytes for 3 bit per every 4-byte word in a 64 byte cache line
    free(values);
    values = NULL;
    
    if(comFlag)
        Word2ByteHuff(request, flag, size, comSize, words, wordPos,huffbit);
    
    
    
    return comFlag;


    

}



bool FRFCFS::HFPCCompress(NVMainRequest *request, uint64_t size, bool flag){
    //GetHuffCode,对64B的cacheline ，16*32bits 每32bits 采样一下，总共采样16个看pattern频率
    // 输出 all zero/8bits符号扩展/16bits符号扩展/32bits/00ab00cd类/重复类/不可压缩类 
    //here outside we still trans size/4 ,so here size*4 means size outside
    //then values[i] is a 32bits thing
    uint64_t * values = convertByte2Word(request, flag, size*4, 8);
    uint64_t i;
    
    uint64_t words[16];
    uint64_t wordPos[16]; //0~8 chars
    uint64_t comSize = 0;
    bool comFlag = false;
    
    if( isZeroPackable( values, size*4 / 8))
    {
        // 000
        // 000 静态全0压缩
        // std::cout<<"64 all zeros" <<std::endl;
        words[0] = 0;
        wordPos[0] = 1;
        comFlag = true;
        comSize = 1;
        free(values);
        values = NULL;
        Word2Byte(request, flag, comSize, comSize, words, wordPos);
        return comFlag;
    }
    free(values);
    values = convertByte2Word(request, flag, size*4, 4);
    int huffbit=0;
    uint64_t Huffreq[8]={0};
    //增加4bits的情况
    for (i = 0; i < size; i++)
    {
        if(values[i]==0){
            Huffreq[0]++;
            continue;
        }
        if(my_abs((int)(values[i])) <= 0xFF){
            Huffreq[1]++;
            continue;
        }
        if(my_abs((int)(values[i])) <= 0xFFFF){
            Huffreq[2]++;
            continue;
        }
        if(((values[i]) & 0xFFFF) == 0 ){
            Huffreq[3]++;
            continue;
        }
        if( my_abs((int)((values[i]) & 0xFFFF)) <= 0xFF
             && my_abs((int)((values[i] >> 16) & 0xFFFF)) <= 0xFF){
            Huffreq[4]++;
            continue;
        }
        if(my_abs((int)(values[i])) <= 0xF){
            Huffreq[7]++;
            continue;
        }
        uint64_t byte0 = (values[i]) & 0xFF;
        uint64_t byte1 = (values[i] >> 8) & 0xFF;
        uint64_t byte2 = (values[i] >> 16) & 0xFF;
        uint64_t byte3 = (values[i] >> 24) & 0xFF;
        if(byte0 == byte1 && byte0 == byte2 && byte0 == byte3){
            Huffreq[5]++;
            continue;
        }
        Huffreq[6]++;
    }
    // for (int i = 0; i < 7; i++)
    // {
    //     std::cout<<"huffmancode["<<i<<"] is"<< Huffreq[i]<<"\t";
    // }
    // printf("\n");
    //h对应新加的4bits
    char name[8] = {'a','b','c','d','e','f','g','h'};
    
    map<char, uint64_t> mapCh;

    for (int i = 0; i < 8; i++)
    {
        mapCh.insert(map<char,uint64_t>::value_type(name[i],Huffreq[i]));
    }
    HuffmanTree huffTree;

    huffTree.Input(mapCh);
    //此时有用的就是veccode里面的。
    // for (int i = 0; i < 7; i++) {
    //     uint64_t temp = strtol(huffTree.veccode[name[i]].c_str(),NULL,10);
    //     std::cout << temp <<"and bits are" <<huffTree.veccode[name[i]].size() <<std::endl;
    // }


    for (i=0;i<size;i++){
        //这里word[i]先只存实际值,用huffbit来存使用的bits数
        //00 -> all zero
        if(values[i] == 0){
            // words[i] = values[i] + 0x0;
            words[i] = values[i] + strtol(huffTree.veccode[name[0]].c_str(),NULL,2);
            // wordPos[i] = 1;
            wordPos[i] = 0;
            comSize += wordPos[i];
            // huffbit +=2;
            huffbit += huffTree.veccode[name[0]].size();
            continue;
        }

        // 01 8bits符号扩展
        if(my_abs((int)(values[i])) <= 0xFF){
            words[i] = my_abs((int)(values[i])) + (strtol(huffTree.veccode[name[1]].c_str(),NULL,2)<<8);
            wordPos[i] = 3;
            comSize += wordPos[i];
            // huffbit +=2;
            huffbit += huffTree.veccode[name[1]].size();
            continue;
        }
        if(my_abs((int)(values[i])) <= 0xF){
            words[i] = my_abs((int)(values[i])) + (strtol(huffTree.veccode[name[7]].c_str(),NULL,2)<<4);
            wordPos[i] = 2;
            comSize += wordPos[i];
            // huffbit +=2;
            huffbit += huffTree.veccode[name[7]].size();
            continue;
        }
        // 110 16bits符号扩展
        if(my_abs((int)(values[i])) <= 0xFFFF){
            // words[i] = my_abs((int)(values[i])) + 0x40000;
            words[i] = my_abs((int)(values[i])) + (strtol(huffTree.veccode[name[2]].c_str(),NULL,2)<<16);
            wordPos[i] = 4;
            comSize += wordPos[i];
            // huffbit +=3;
            huffbit += huffTree.veccode[name[2]].size();
            continue;
        }
        //100  32bits
        if(((values[i]) & 0xFFFF) == 0 ){
            // words[i] = (values[i] >> 16) + 0x40000;
            words[i] = (values[i] >> 16) + (strtol(huffTree.veccode[name[3]].c_str(),NULL,2)<<16);
            wordPos[i] = 4;
            comSize += wordPos[i];
            // huffbit +=3;
            huffbit += huffTree.veccode[name[3]].size();
            continue;
        }
        //101 00ab00cd类
        if( my_abs((int)((values[i]) & 0xFFFF)) <= 0xFF
             && my_abs((int)((values[i] >> 16) & 0xFFFF)) <= 0xFF){
            words[i] = my_abs((int)((values[i] >> 8))) + my_abs((int)((values[i]) & 0xFFFF)) + (strtol(huffTree.veccode[name[4]].c_str(),NULL,2)<<32);
            wordPos[i] = 4;
            comSize += wordPos[i];
            // huffbit +=3;
            huffbit += huffTree.veccode[name[4]].size();
            continue;
        }
        //1110 重复类
        uint64_t byte0 = (values[i]) & 0xFF;
        uint64_t byte1 = (values[i] >> 8) & 0xFF;
        uint64_t byte2 = (values[i] >> 16) & 0xFF;
        uint64_t byte3 = (values[i] >> 24) & 0xFF;
        if(byte0 == byte1 && byte0 == byte2 && byte0 == byte3){
            words[i] = byte0 +  (strtol(huffTree.veccode[name[5]].c_str(),NULL,2)<<8);
            wordPos[i] = 2;
            comSize += wordPos[i];
            // huffbit +=4;
            huffbit += huffTree.veccode[name[5]].size();
            continue;
        }
        //1111
        words[i] = values[i];
        wordPos[i] = 8;
        //这里huffbit不用加因为有一bit是是否压缩的
        comSize += wordPos[i];
    }
    //这里我们把prefix加上去
    // std::cout << "now huffbit is"<<huffbit <<std::endl;
    comSize += huffbit/4;
    if (comSize==0)
        comSize=1;
    if(comSize % 2 == 1)
        comSize++;
    comSize /= 2;
    if(comSize < (size*4))
    {
        comFlag = true;
    }
    
    //6 bytes for 3 bit per every 4-byte word in a 64 byte cache line
    free(values);
    values = NULL;
    
    if(comFlag)
        Word2ByteHuff(request, flag, size, comSize, words, wordPos,huffbit);
    
    
    
    return comFlag;


    

}


bool FRFCFS::FPCCompress(NVMainRequest *request, uint64_t size, bool flag ){
    uint64_t * values = convertByte2Word(request, flag, size*4, 4);
    uint64_t i;
    uint64_t words[16];
    uint64_t wordPos[16]; //0~8 chars
    uint64_t comSize = 0;
    bool comFlag = false;
    
    for (i = 0; i < size; i++) {
     
        // 000
        if(values[i] == 0){
            words[i] = values[i] + 0x0;
            wordPos[i] = 1;
            comSize += wordPos[i];
            continue;
        }
        // 001
        if(my_abs((int)(values[i])) <= 0xFF){
            words[i] = my_abs((int)(values[i])) + 0x100;
            wordPos[i] = 3;
            comSize += wordPos[i];
            continue;
        }
        // 011
        if(my_abs((int)(values[i])) <= 0xFFFF){
            words[i] = my_abs((int)(values[i])) + 0x30000;
            wordPos[i] = 5;
            comSize += wordPos[i];
            continue;
        }
        //100  
        if(((values[i]) & 0xFFFF) == 0 ){
            words[i] = (values[i] >> 16) + 0x40000;
            wordPos[i] = 5;
            comSize += wordPos[i];
            continue;
        }
        //101
        if( my_abs((int)((values[i]) & 0xFFFF)) <= 0xFF
             && my_abs((int)((values[i] >> 16) & 0xFFFF)) <= 0xFF){
            words[i] = my_abs((int)((values[i] >> 8))) + my_abs((int)((values[i]) & 0xFFFF)) + 0x50000;
            wordPos[i] = 5;
            comSize += wordPos[i];
            continue;
        }
        //110
        uint64_t byte0 = (values[i]) & 0xFF;
        uint64_t byte1 = (values[i] >> 8) & 0xFF;
        uint64_t byte2 = (values[i] >> 16) & 0xFF;
        uint64_t byte3 = (values[i] >> 24) & 0xFF;
        if(byte0 == byte1 && byte0 == byte2 && byte0 == byte3){
            words[i] = byte0 + 0x600;
            wordPos[i] = 3;
            comSize += wordPos[i];
            continue;
        }
        //111
        words[i] = values[i];
        wordPos[i] = 8;
        comSize += wordPos[i];
    }
    if(comSize % 2 == 1)
        comSize++;
    comSize /= 2;
    if(comSize < (size*4))
    {
        comFlag = true;
    }
    
    //6 bytes for 3 bit per every 4-byte word in a 64 byte cache line
    free(values);
    values = NULL;
    
    if(comFlag)
        Word2Byte(request, flag, size, comSize, words, wordPos);
    
    
    
    return comFlag;
        
}

bool FRFCFS::isZeroPackable ( uint64_t * values, uint64_t size){
    bool nonZero = false;
    uint64_t i;
    for (i = 0; i < size; i++) {
        if( values[i] != 0){
            nonZero = true;
            break;
        }
    }
    return !nonZero;
}

///
/// Check if the cache line consists of only same values
///
bool FRFCFS::isSameValuePackable ( uint64_t * values, uint64_t size){
    bool notSame = false;
    uint64_t i;
    for (i = 0; i < size; i++) {
        if( values[0] != values[i]){
            notSame = true;
            break;
        }
    }
    return !notSame;
}

uint64_t FRFCFS::multBaseCompression ( uint64_t * values, uint64_t size, uint64_t blimit, uint64_t bsize, uint64_t *currWords, uint64_t *currWordPos, uint64_t &pos)
{
    uint64_t limit = 0;
    uint64_t BASES = 2;
    //define the appropriate size for the mask
    //size = 8,blimit = 1 , bsize = 8
    switch(blimit){
        case 1:
            limit = 0xFF;
            break;
        case 2:
            limit = 0xFFFF;
            break;
        case 4:
            limit = 0xFFFFFFFF;
            break;
        default:
            break;

    }
    uint64_t mbases[2]={0};
    uint64_t baseCount = 1;
    mbases[0] = 0;
    uint64_t i,j;
    for (i = 0; i < size; i++) {
        //long long int 64bits
        //只要有一个values的值和base的区别大，就设置为第二个基
        if( my_llabs((long long int)(mbases[0] -  values[i])) > limit ){
            // add new base
            mbases[1] = values[i];
            baseCount++;
            break;
        }
            
    }
    // find how many elements can be compressed with mbases
    uint64_t compCount = 0;
    //这里的bsize是BASE的size，如果大于4,nums为2
    //bsize = 8 ,nums=2 curwords[0]=
    uint64_t nums = (bsize > 4)?2:1;
    //对于每个基来说
    for (pos = 0; pos < baseCount; pos++)
    {
        for(i = 0; i < nums; i++)
        //对于每个基的nums来说
        {
            currWords[pos*nums + i] = (values[pos] >> (32*(1-i))) & 0xFFFFFFFF;
            currWordPos[pos*nums + i] = (bsize*2 > 8)?8:(bsize*2);
            //例如这里currwords就是values[0/1]的高32，低32位
        }
    }
    pos = pos * nums;
    //pos = 2*pos
    for (i = 0; i < size; i++) {
        for(j = 0; j <  baseCount; j++){
            //对于每个基而言，如果相减小于limit当前字为delta，且占位为blimit的两倍，blimit=1 相当于limit = 0xff
            if( my_llabs((long long int)(mbases[j] -  values[i])) <= limit ){
                //limit * 2
                currWords[pos] = my_llabs((long long int)(mbases[j] -  values[i])) & limit;
                currWordPos[pos++] = blimit * 2;
                //blimit = 1 说明用了0xff共8位去压缩， 8/4 = 2 = 1*2
                compCount++;
                break;
            }
        }
    }
    //return compressed size
    //总共有size个，压缩了compCount个，压缩的大小为blimit，未压缩的大小为bsize
    uint64_t mCompSize = blimit * compCount + bsize * BASES + (size - compCount) * bsize;
    if(compCount < size)
        return size * bsize;
    return mCompSize;
}

bool FRFCFS::BDICompress (NVMainRequest *request, uint64_t _blockSize, bool flag )
{
    //blocksize就是req的data的size，一般是64
    uint64_t * values = convertByte2Word(request, flag, _blockSize, 8);
    //把64个指向uint8的转换成8个word,一个word64bits
    //也就是论文中说的8*8bytes
    uint64_t bestCSize = _blockSize;
    uint64_t currCSize = _blockSize;
    uint64_t i, pos, bestPos;
    uint64_t words[35];
    uint64_t wordPos[35]; //0~8 chars
    
    uint64_t currWords[35];
    uint64_t currWordPos[35]; //0~8 chars
    bool comFlag = false;
    bestPos = 16;
    //看一条cacheline内的差别，比较value[0]~value[8]
    if( isSameValuePackable( values, _blockSize / 8))
    {
        currCSize = 8;
    }
    if(bestCSize > currCSize)
    {
        //这里我也不太清楚到底是哪种，大概是全0或者全相同中的一种
        //bestPos就=2，word[1]为value高32bits，word[2]低32bits
        bestCSize = currCSize;
        bestPos = bestCSize / 4;
        //假设这个value 1-8是一样的, bestPos = 8/4 =2
        words[0] = 0x0;
        wordPos[0] = 1;
        for(i = 0; i < bestPos; i++)
        {
            words[i+1] = (values[i/2] >> (32*(1-i%2))) & 0xFFFFFFFF;
            wordPos[i+1] = 8;
        }
        //words[1] = value[0]高32位
        //words[2] = value[0]低32位
        bestPos++;
        //然后bestPos=3
    }
    //传入values,8,1,8,....
    //总共有size = 8个，blimit=1,bsize=8为不压缩的大小
    currCSize = multBaseCompression( values, _blockSize / 8, 1, 8, currWords, currWordPos, pos);
    if(bestCSize > currCSize)
    {
        bestCSize = currCSize;
        bestPos = pos;
        words[0] = 0x1;
        wordPos[0] = 1;
        for(i = 0; i < bestPos; i++)
        {
            words[i+1] = currWords[i];
            wordPos[i+1] = currWordPos[i];
        }
        bestPos++;
    }
    currCSize = multBaseCompression( values, _blockSize / 8, 2, 8, currWords, currWordPos, pos);
    if(bestCSize > currCSize)
    {
        bestCSize = currCSize;
        bestPos = pos;
        words[0] = 0x2;
        wordPos[0] = 1;
        for(i = 0; i < bestPos; i++)
        {
            words[i+1] = currWords[i];
            wordPos[i+1] = currWordPos[i];
        }
        bestPos++;
    }
    currCSize = multBaseCompression( values, _blockSize / 8, 4, 8, currWords, currWordPos, pos);
    if(bestCSize > currCSize)
    {
        bestCSize = currCSize;
        bestPos = pos;
        words[0] = 0x3;
        wordPos[0] = 1;
        for(i = 0; i < bestPos; i++)
        {
            words[i+1] = currWords[i];
            wordPos[i+1] = currWordPos[i];
        }
        bestPos++;
    }
    free(values);
    values = convertByte2Word(request, flag, _blockSize, 4);
    if( isSameValuePackable( values, _blockSize / 4))
    {
        currCSize = 4;
    }
    if(bestCSize > currCSize)
    {
        bestCSize = currCSize;
        bestPos = bestCSize / 4;
        words[0] = 0x4;
        wordPos[0] = 1;
        for(i = 0; i < bestPos; i++)
        {
            words[i+1] = (values[i/2] >> (32*(1-i%2))) & 0xFFFFFFFF;
            wordPos[i+1] = 8;
        }
        bestPos++;
    }
    currCSize = multBaseCompression( values, _blockSize / 4, 1, 4, currWords, currWordPos, pos);
    if(bestCSize > currCSize)
    {
        bestCSize = currCSize;
        bestPos = pos;
        words[0] = 0x5;
        wordPos[0] = 1;
        for(i = 0; i < bestPos; i++)
        {
            words[i+1] = currWords[i];
            wordPos[i+1] = currWordPos[i];
        }
        bestPos++;
    }
    currCSize = multBaseCompression( values, _blockSize / 4, 2, 4, currWords, currWordPos, pos);
    if(bestCSize > currCSize)
    {
        bestCSize = currCSize;
        bestPos = pos;
        words[0] = 0x6;
        wordPos[0] = 1;
        for(i = 0; i < bestPos; i++)
        {
            words[i+1] = currWords[i];
            wordPos[i+1] = currWordPos[i];
        }
        bestPos++;
    }
    free(values);
    values = convertByte2Word(request, flag, _blockSize, 2);
    currCSize = multBaseCompression( values, _blockSize / 2, 1, 2, currWords, currWordPos, pos);
    if(bestCSize > currCSize)
    {
        bestCSize = currCSize;
        bestPos = pos;
        words[0] = 0x7;
        wordPos[0] = 1;
        for(i = 0; i < bestPos; i++)
        {
            words[i+1] = currWords[i];
            wordPos[i+1] = currWordPos[i];
        }
        bestPos++;
    }
    free(values);
    values = NULL;
    
    if(bestCSize < _blockSize)
    {
        comFlag = true;
        Word2Byte(request, flag, bestPos, bestCSize, words, wordPos);
    }
    
    return comFlag;

}

bool FRFCFS::DFPCCompress(NVMainRequest *request, uint64_t _blockSize )
{
    if(mem_writes < granularities)
    //如果写次数小于粒度
	{
        //所以这里为什么要除以4呢
        FPCIdentify(request, _blockSize / 4);
        BDIIdentify(request, _blockSize);
        Sample(request, _blockSize);
        StaticCompress(request, _blockSize/4, false);
		return StaticCompress(request, _blockSize/4, true);
	}else
	{
        if(sample_flag)
            ExtractPattern();
        DynamicCompress(request, _blockSize, false);
		return DynamicCompress(request, _blockSize, true);
	}
}

bool FRFCFS::HDFPCCompress(NVMainRequest *request, uint64_t _blockSize )
{
    if(mem_writes < granularities)
    //如果写次数小于粒度
	{
        //所以这里为什么要除以4呢
        FPCIdentify(request, _blockSize / 4);
        BDIIdentify(request, _blockSize);
        Sample(request, _blockSize);
        // StaticCompress(request, _blockSize/4, false);

        HFPCCompress(request,_blockSize/4,false);
		// return StaticCompress(request, _blockSize/4, true);
        return HFPCCompress(request,_blockSize/4,true);
	}else
	{
        if (buildtree==false){
            printf("我们终于新建了一棵树\n");
            buildtree = BuildHuffTree(huffTree1,mapCh1);
            //这里我们新建一棵树
        }
        if(sample_flag)
            ExtractPattern();
        // printf("进入DynamicHuffCompress\n");
        DynamicHuffCompress(request, _blockSize, false);
		return DynamicHuffCompress(request, _blockSize, true);
	}
}

bool FRFCFS::BuildHuffTree(HuffmanTree &huffTree1,map<char, uint64_t> &mapCh1)
{
    char name[8] = {'a','b','c','d','e','f','g','h'};
    
    // map<char, uint64_t> mapCh1;

    for (int i = 0; i < 8; i++)
    {   
        std::cout << "pattern"<<pattern_ana[i] <<std::endl;
        mapCh1.insert(map<char,uint64_t>::value_type(name[i],pattern_ana[i]));
    }
    // HuffmanTree huffTree1;

    huffTree1.Input(mapCh1);

    return true;
}



//这里还是求了fpcsize的，如果不求的话或许可以降低延迟？
bool FRFCFS::DynamicHuffCompress(NVMainRequest *request, uint64_t size, bool flag  )
{
    
    uint64_t FPC_pattern_size=0, BDI_pattern_size=0,HFPC_pattern_size=0;
    // printf("进入了DynamicHuffCompress\n");
    
    FPC_pattern_size = DynamicFPCCompress(request, size/4, flag);
    //pure就是那种几千万次做的hufftree
    HFPC_pattern_size = PUREHFPCCompress(request,size/4,flag);
    if(FPC_pattern_size<HFPC_pattern_size)
    {
        printf("fpc is better than hfpc?");
    }
    BDI_pattern_size = DynamicBDICompress(request, size, flag);
    if(HFPC_pattern_size < BDI_pattern_size)
    {
        PUREHFPCCompress(request,size/4,flag);
    }
    if(flag)
        return request->data.IsCompressed();
    else
        return request->oldData.IsCompressed();
    
}

bool FRFCFS::DynamicCompress(NVMainRequest *request, uint64_t size, bool flag  )
{
    
    uint64_t FPC_pattern_size=0, BDI_pattern_size=0;
    
    FPC_pattern_size = DynamicFPCCompress(request, size/4, flag);
    BDI_pattern_size = DynamicBDICompress(request, size, flag);
    if(FPC_pattern_size < BDI_pattern_size)
    {
        DynamicFPCCompress(request, size/4, flag);
    }
    if(flag)
        return request->data.IsCompressed();
    else
        return request->oldData.IsCompressed();
    
}

uint64_t FRFCFS::DynamicFPCCompress(NVMainRequest *request, uint64_t size, bool flag )
{
    uint64_t * values = convertByte2Word(request, flag, size*4, 8);
    //这里传进来的就是size/4,乘了4以后又变成了 values指向req里面指向uint8的类型的个数
    //所以这里就变成了
    uint64_t i;
    
    uint64_t words[16];
    uint64_t wordPos[16]; //0~8 chars
    uint64_t comSize = 0;
    bool comFlag = false;
    bool dynamicFlag = false;
    
    if( isZeroPackable( values, size*4 / 8))
    {
        // 000
        //这里是把每个req划分成64bits，若都为0
        words[0] = 0;
        wordPos[0] = 1;
        comFlag = true;
        comSize = 1;
        free(values);
        values = NULL;
        //这里的word2byte是没办法继续压缩的
        Word2Byte(request, flag, comSize, comSize, words, wordPos);
        //这里req就直接返回了，所以req改了也没关系
        
        return comSize;
    }
    free(values);
    //重新，这时候step=4
    values = convertByte2Word(request, flag, size*4, 4);
    //下面这个size实际上就是实际的32位的个数，就是有多少个values，因为 size是/4过的
    for (i = 0; i < size; i++) {
        // 001
        if(values[i] == 0){
            words[i] = values[i] + 0x1;
            wordPos[i] = 1;
            comSize += wordPos[i];
            continue;
        }
        dynamicFlag = true;
        //默认mask_pos=0
        //这里的两个函数要等提取模式看了再回来看到底有什么用
        for(int j = 0; j < mask_pos && dynamicFlag; j++)
        {
            int compressible_char = 8 - compressibleChars[j];
            if(compressible_char < 4 && ((my_abs((int)(values[i])) & masks[j]) == 0))
            {
                words[i] = ((j+4)<<(compressible_char*4));
                uint64_t mask = masks[j];
                uint64_t word = my_abs((int)(values[i]));
                for(int k = 0; k < compressible_char; k++)
                {
                    while((mask & 0xF) != 0)
                    {
                        mask = mask >> 4;
                        word = word >> 4;
                    }
                    words[i] = words[i] | ((word & 0xF) << (k * 4));
                }
                //words[i] = my_abs((int)(values[i])) + ((j+4)<<(compressible_char*4));

                wordPos[i] = 1 + compressible_char;
                comSize += wordPos[i];
                dynamicFlag = false;
                break;
            }
        }
        if(!dynamicFlag)
            continue;
        // 011
        //低16位的
        if(my_abs((int)(values[i])) <= 0xFFFF){
            words[i] = my_abs((int)(values[i])) + 0x30000;
            wordPos[i] = 5;
            comSize += wordPos[i];
            continue;
        }
        //100  
        // 高16位的
        if(((values[i]) & 0xFFFF) == 0 ){
            words[i] = (values[i] >> 16) + 0x40000;
            wordPos[i] = 5;
            comSize += wordPos[i];
            continue;
        }
        for(int j = 0; j < mask_pos && dynamicFlag; j++)
        {
            int compressible_char = 8 - compressibleChars[j];
            if(compressible_char >= 4 && ((my_abs((int)(values[i])) & masks[j]) == 0))
            {
                
                words[i] = ((j+4)<<(compressible_char*4));
                uint64_t mask = masks[j];
                uint64_t word = my_abs((int)(values[i]));
                for(int k = 0; k < compressible_char; k++)
                {
                    while((mask & 0xF) != 0)
                    {
                        mask = mask >> 4;
                        word = word >> 4;
                    }
                    words[i] = words[i] | ((word & 0xF) << (k * 4));
                }
                //words[i] = my_abs((int)(values[i])) + ((j+4)<<(compressible_char*4));
                wordPos[i] = 1 + compressible_char;
                //这样真的有压缩过吗？
                comSize += wordPos[i];
                dynamicFlag = false;
                break;
            }
        }
        if(!dynamicFlag)
            continue;
        
        //110
        //这个是有1+BDICOUNT个的bool型数组,default = false
        //110代表重复aa/aa/aa/aa
        if(special_pattern_flag[0])
        {
            uint64_t byte0 = (values[i]) & 0xFF;
            uint64_t byte1 = (values[i] >> 8) & 0xFF;
            uint64_t byte2 = (values[i] >> 16) & 0xFF;
            uint64_t byte3 = (values[i] >> 24) & 0xFF;
            if(byte0 == byte1 && byte0 == byte2 && byte0 == byte3){
                words[i] = byte0 + 0x600;
            wordPos[i] = 3;
            comSize += wordPos[i];
                continue;
            }
        }
        //上面这些都不行就是压缩不了
        words[i] = values[i];
        wordPos[i] = 8;
        comSize += wordPos[i];
    }
    //因为一直都是4bits 4bits的来的
    if(comSize % 2 == 1)
        comSize++;
    comSize /= 2;
    if(comSize < (size*4))
    {
        comFlag = true;
    }
    free(values);
    values = NULL;
    if(comFlag)
    {
        Word2Byte(request, flag, size, comSize, words, wordPos);
        if(flag)
            comSize = request->data.GetComSize();
    }
    
    
    //6 bytes for 3 bit per every 4-byte word in a 64 byte cache line
    
    return comSize;
}

uint64_t FRFCFS::DynamicBDICompress(NVMainRequest *request, uint64_t _blockSize, bool flag )
{
    uint64_t * values = convertByte2Word(request, flag, _blockSize, 8);
    uint64_t bestCSize = _blockSize;
    uint64_t currCSize = _blockSize;
    uint64_t i, pos, bestPos;
    uint64_t words[34];
    uint64_t wordPos[34]; //0~8 chars
    
    uint64_t currWords[34];
    uint64_t currWordPos[34]; //0~8 chars
    bestPos = 0;
    //这里的special_pattern_flag用处暂时不清楚
    //好像之所以是BDI+1,因为[0]用在FPC上了
    if(special_pattern_flag[1] || special_pattern_flag[2] || special_pattern_flag[3] || special_pattern_flag[4])
    {
        if( isSameValuePackable( values, _blockSize / 8))
        {
            currCSize = 8;
        }
        if(bestCSize > currCSize)
        {
            bestCSize = currCSize;
            bestPos = bestCSize / 4;
            for(i = 0; i < bestPos; i++)
            {
                words[i] = (values[i/2] >> (32*(1-i%2))) & 0xFFFFFFFF;
                wordPos[i] = 8;
                //这里的8个实际上是32bits
            }
        }
        currCSize = multBaseCompression( values, _blockSize / 8, 1, 8, currWords, currWordPos, pos);
        if(bestCSize > currCSize)
        {
            bestCSize = currCSize;
            bestPos = pos;
            for(i = 0; i < bestPos; i++)
            {
                words[i] = currWords[i];
                wordPos[i] = currWordPos[i];
            }
        }
        currCSize = multBaseCompression( values, _blockSize / 8, 2, 8, currWords, currWordPos, pos);
        if(bestCSize > currCSize)
        {
            bestCSize = currCSize;
            bestPos = pos;
            for(i = 0; i < bestPos; i++)
            {
                words[i] = currWords[i];
                wordPos[i] = currWordPos[i];
            }
        }
        currCSize = multBaseCompression( values, _blockSize / 8, 4, 8, currWords, currWordPos, pos);
        if(bestCSize > currCSize)
        {
            bestCSize = currCSize;
            bestPos = pos;
            for(i = 0; i < bestPos; i++)
            {
                words[i] = currWords[i];
                wordPos[i] = currWordPos[i];
            }
        }
    }
    if(special_pattern_flag[5] || special_pattern_flag[6] || special_pattern_flag[7])
    {
        free(values);
        values = convertByte2Word(request, flag, _blockSize, 4);
        if( isSameValuePackable( values, _blockSize / 4))
        {
            currCSize = 4;
        }
        if(bestCSize > currCSize)
        {
            bestCSize = currCSize;
            bestPos = bestCSize / 4;
            for(i = 0; i < bestPos; i++)
            {
                words[i] = (values[i/2] >> (32*(1-i%2))) & 0xFFFFFFFF;
                wordPos[i] = 8;
            }
        }
        
        currCSize = multBaseCompression( values, _blockSize / 4, 1, 4, currWords, currWordPos, pos);
        if(bestCSize > currCSize)
        {
            bestCSize = currCSize;
            bestPos = pos;
            for(i = 0; i < bestPos; i++)
            {
                words[i] = currWords[i];
                wordPos[i] = currWordPos[i];
            }
        }
        currCSize = multBaseCompression( values, _blockSize / 4, 2, 4, currWords, currWordPos, pos);
        if(bestCSize > currCSize)
        {
            bestCSize = currCSize;
            bestPos = pos;
            for(i = 0; i < bestPos; i++)
            {
                words[i] = currWords[i];
                wordPos[i] = currWordPos[i];
            }
        }
    }
    if(special_pattern_flag[8])
    {
        free(values);
        values = convertByte2Word(request, flag, _blockSize, 2);
        currCSize = multBaseCompression( values, _blockSize / 2, 1, 2, currWords, currWordPos, pos);
        if(bestCSize > currCSize)
        {
            bestCSize = currCSize;
            bestPos = pos;
            for(i = 0; i < bestPos; i++)
            {
                words[i] = currWords[i];
                wordPos[i] = currWordPos[i];
            }
        }
    }
    free(values);
    values = NULL;
    if(bestCSize < _blockSize)
    {
        //word2byte第一步就是 把comsize设置成bestcsize，然后如果可以压缩就改小comsize
        Word2Byte(request, flag, bestPos, bestCSize, words, wordPos);
        if(flag)
            bestCSize = request->data.GetComSize();
        /*else
            bestCSize = request->oldData.GetComSize();*/
    }
    
    
    return bestCSize;
}

bool FRFCFS::StaticCompress(NVMainRequest *request, uint64_t size, bool flag )
{
    uint64_t * values = convertByte2Word(request, flag, size*4, 8);
    uint64_t i;
    
    uint64_t words[16];
    uint64_t wordPos[16]; //0~8 chars
    uint64_t comSize = 0;
    bool comFlag = false;
    
    if( isZeroPackable( values, size*4 / 8))
    {
        // 000
        // 000 静态全0压缩
        words[0] = 0;
        wordPos[0] = 1;
        comFlag = true;
        comSize = 1;
        free(values);
        values = NULL;
        Word2Byte(request, flag, comSize, comSize, words, wordPos);
        return comFlag;
    }
    free(values);
    values = convertByte2Word(request, flag, size*4, 4);
    for (i = 0; i < size; i++) {
     
        // 001
        if(values[i] == 0){
            words[i] = values[i] + 0x1;
            wordPos[i] = 1;
            comSize += wordPos[i];
            continue;
        }
        // 010
        if(my_abs((int)(values[i])) <= 0xFFFF){
            words[i] = my_abs((int)(values[i])) + 0x20000;
            wordPos[i] = 5;
            comSize += wordPos[i];
            continue;
        }
        //011  
        if(((values[i]) & 0xFFFF) == 0 ){
            words[i] = (values[i] >> 16) + 0x30000;
            wordPos[i] = 5;
            comSize += wordPos[i];
            continue;
        }
        //uncompressible
        words[i] = values[i];
        wordPos[i] = 8;
        comSize += wordPos[i];
    }
    if(comSize % 2 == 1)
        comSize++;
    comSize /= 2;
    if(comSize < (size*4))
    {
        comFlag = true;
    }
    if(comFlag)
        Word2Byte(request, flag, size, comSize, words, wordPos);
    free(values);
    values = NULL;
    return comFlag;
}

uint64_t FRFCFS::FPCIdentify(NVMainRequest *request, uint64_t size){
    uint64_t * values = convertByte2Word(request, true, size*4, 4);
    uint64_t i;
    for (i = 0; i < size; i++) {
        if(values[i] == 0){
            continue;
        }
        if(my_abs((int)(values[i])) <= 0xFF){
			FPCCounter[0]++;
            continue;
        }
        // 011
        if(my_abs((int)(values[i])) <= 0xFFFF){
			
            continue;
        }
        //100  
        if(((values[i]) & 0xFFFF) == 0 ){
			
            continue;
        }
        //101
        //00ab00cc有补码
        if( my_abs((int)((values[i]) & 0xFFFF)) <= 0xFF
             && my_abs((int)((values[i] >> 16) & 0xFFFF)) <= 0xFF){
            FPCCounter[1]++;
			
            continue;
        }
        //110
        uint64_t byte0 = (values[i]) & 0xFF;
        uint64_t byte1 = (values[i] >> 8) & 0xFF;
        uint64_t byte2 = (values[i] >> 16) & 0xFF;
        uint64_t byte3 = (values[i] >> 24) & 0xFF;
        if(byte0 == byte1 && byte0 == byte2 && byte0 == byte3){
			FPCCounter[2]++;
			
            continue;
        }
        
    }
    free(values);
    values = NULL;
    return 1;
        
}

uint64_t FRFCFS::BDIIdentify (NVMainRequest *request, uint64_t _blockSize)
{
 
    uint64_t * values = convertByte2Word(request, true, _blockSize, 8);
    uint64_t bestCSize = _blockSize;
    uint64_t currCSize = _blockSize;
	bool isBest[8];
	uint64_t currWords[34];
    uint64_t currWordPos[34]; //0~8 chars
    uint64_t pos;
	for(int i=0;i<8;i++)
	{
		isBest[i] = false;
	}
    if( isSameValuePackable( values, _blockSize / 8))
    {
        currCSize = 8;
		isBest[0] = true;
    }
    if(bestCSize > currCSize)
    {
        bestCSize = currCSize;
    }
    currCSize = multBaseCompression( values, _blockSize / 8, 1, 8, currWords, currWordPos, pos);
    if(bestCSize > currCSize)
    {
        bestCSize = currCSize;
		isBest[1] = true;
    }
    currCSize = multBaseCompression( values, _blockSize / 8, 2, 8, currWords, currWordPos, pos);
    if(bestCSize > currCSize)
    {
        bestCSize = currCSize;
		isBest[2] = true;
    }
    currCSize = multBaseCompression( values, _blockSize / 8, 4, 8, currWords, currWordPos, pos);
    if(bestCSize > currCSize)
    {
        bestCSize = currCSize;
		isBest[3] = true;
    }
    free(values);
    values = convertByte2Word(request, true, _blockSize, 4);
    if( isSameValuePackable( values, _blockSize / 4))
    {
        currCSize = 4;
    }
    if(bestCSize > currCSize)
    {
        bestCSize = currCSize;
        isBest[4] = true;
    }
    currCSize = multBaseCompression( values, _blockSize / 4, 1, 4, currWords, currWordPos, pos);
    if(bestCSize > currCSize)
    {
        bestCSize = currCSize;
        isBest[5] = true;
    }
    currCSize = multBaseCompression( values, _blockSize / 4, 2, 4, currWords, currWordPos, pos);
    if(bestCSize > currCSize)
    {
        bestCSize = currCSize;
        isBest[6] = true;
    }
    free(values);
    values = convertByte2Word(request, true, _blockSize, 2);
    currCSize = multBaseCompression( values, _blockSize / 2, 1, 2, currWords, currWordPos, pos);
    if(bestCSize > currCSize)
    {
        bestCSize = currCSize;
        isBest[7] = true;
    }
    free(values);
    values = NULL;
	for(int i=7;i>=0;i--)
	{
		if(isBest[i])
		{
			BDICounter[i] += _blockSize - bestCSize;
            BDIpatterncounter++;
            break;
		}
	}
    return 1;

}




//这个sampler是对64B的cacheline ，16*32bits 每4bits 采样一下，总共采样128个看哪些为0
uint64_t FRFCFS::Sample(NVMainRequest *request, uint64_t _blockSize)
{
	uint64_t num = 0;
	uint64_t * values = convertByte2Word(request, true, _blockSize, 4);
	for(int i = 0; i< 16;i++)
	{
		num = my_llabs((long long int)values[i]);
		for(int j=0; j<8;j++)
		{
			if((num & 0xF) == 0)
			{
				SampleCounter[ 8*i+j ]++;
			}
			num = (num >> 4);
		}
	}
    free(values);
    values = NULL;
	return 1;
}

uint64_t FRFCFS::ExtractPattern()
{
    int i, j, pos;
    bool flag;
    uint64_t lower_bound, upper_bound, threshold;//, word_size;
    int word_count = SAMPLECOUNT/DYNAMICWORDSIZE; // 32-bit word
    //128 / 8 = 16 因为一个字是32位就是 8个4bits
    int total_pattern_count = FPCCOUNT+BDICOUNT+word_count;
    
    uint8_t SamplePatterns[SAMPLECOUNT];
    pattern_num = word_count / 2;
    //这个除以二是什么鬼
    uint8_t patterns_temp[word_count];
    int compressible_chars[word_count];
    uint8_t extracted_patterns[word_count];

    uint64_t DynamicPatterns[total_pattern_count];
    uint64_t CompressBytes[total_pattern_count];
    
    
	for(i = 0; i < FPCCOUNT; i++)
    {
        DynamicPatterns[i] = i;
        if(i == 1)
            //CompressBytes[i] = FPCCounter[i] - FPCCounter[i]/2 - FPCCounter[i] * 3 / 32;
            CompressBytes[i] = FPCCounter[i]*4 - FPCCounter[i]*2 - FPCCounter[i] * 3  / 8;
        //(FPCCounter[i] * 4 - FPCCounter[i] * 2 - FPCCounter[i] * 3  / 8) / 4;
        else
            //CompressBytes[i] = FPCCounter[i] - FPCCounter[i]/4 - FPCCounter[i] * 3 / 32;
            CompressBytes[i] = FPCCounter[i]*4 - FPCCounter[i] - FPCCounter[i] * 3 / 8;
        //FPCCompressBytes[i] = FPCCounter[i] * 22;
        //(FPCCounter[i] + 3 * FPCCounter[i] / 8) * 64 / 4;
        // std::cout<<"FPCCompressBytes["<<i<<"]: "<<CompressBytes[i]<<std::endl;
        FPCCounter[i] = 0;
        //提取完了就置零
    }
    
    for(i = 0; i < BDICOUNT; i++)
    {
        CompressBytes[i + FPCCOUNT] = 0;
    }
    
    for(i = 0; i < BDICOUNT; i++)
    {
        DynamicPatterns[i + FPCCOUNT] = i + FPCCOUNT;
        //CompressBytes[i + FPCCOUNT] = (BDICounter[i] - 3 * BDIpatterncounter / 8)/64;
        if(i < 4)
            CompressBytes[0 + FPCCOUNT] += BDICounter[i];
            //这里是前4个除以8的
        else if(i < BDICOUNT - 1)
            CompressBytes[4 + FPCCOUNT] += BDICounter[i];
        else
            CompressBytes[i + FPCCOUNT] = BDICounter[i];
        //std::cout<<"BDICompressBytes["<<i<<"]: "<<CompressBytes[i+FPCCOUNT]<<std::endl;
        BDICounter[i] = 0;
    }
    //    BDICounter[i] = BDICounter[i];
    //一直在思考前两步到底在干嘛，Dynamicpatterns[i] -> i ， 但是compressbytes 只有012，3，7，8有值
    //fpc conter里面存的是32bits符合pattren的个数，但是bdicounter里面村的是压缩节省的size（bits/byte暂时没看）

	
    //sample
    lower_bound = upper_bound = SampleCounter[0];
    
    for(i = 1; i < SAMPLECOUNT; i++)
    {
        if(SampleCounter[i] < lower_bound)
            lower_bound = SampleCounter[i];
        else if(SampleCounter[i] > upper_bound)
            upper_bound = SampleCounter[i];
    }
    ;
    threshold = lower_bound + (upper_bound - lower_bound) * threshold_factor;
    printf("%ld\n", threshold);
    for(i = 0; i < SAMPLECOUNT; i++)
    {
        if(SampleCounter[i] < threshold)
            SamplePatterns[i] = 1;
        else
            SamplePatterns[i] = 0;
        printf("%d\t", SamplePatterns[i]);
    }
    printf("\n");
    //这里根据128个采样的0的数据得到了samplepatterns的值
    //determine word_size
    
    for(i = 0; i < word_count ; i++)
    {
        patterns_temp[i] = 0;
        for(j = 0; j < DYNAMICWORDSIZE; j++)
        {
            patterns_temp[i] = patterns_temp[i] + (SamplePatterns[ DYNAMICWORDSIZE*i+j ] << (7-j));
        }
    }
    for(pattern_num = word_count / 2; pattern_num < word_count && pattern_num >= 1;pattern_num /= 2)
    {
        flag = true;
        for(i = 0; i < pattern_num; i++)
        {
            if(patterns_temp[i] != patterns_temp[i + pattern_num])
            {
                pattern_num *= 2; //restore
                flag = false;
                break;
            }
        }
        if(!flag)
            break;
    }
    if(pattern_num == 0)
        pattern_num = 1;
    printf("pattern_num: %d\n", pattern_num);
    //word_size = pattern_num * 4;//word_size Bytes
    
    //dynamic patterns: extracted_patterns[0~(pos-1)]
    for(pos = 0, i = 0; i < pattern_num; i++)
    {
        uint32_t compressible_char = 0;
        uint64_t min_compression_counter = mem_writes;
        printf("patterns_temp[%d]: %d\t", i, patterns_temp[i]);
        for(j = 0; j < DYNAMICWORDSIZE; j++)
        {
            uint8_t compression_tag = (patterns_temp[i] >> (7 - j)) & 0x1;
            if(compression_tag == 0)
            {
                compressible_char++;
                if(SampleCounter[ DYNAMICWORDSIZE*i+j ] < min_compression_counter)
                    min_compression_counter = SampleCounter[ DYNAMICWORDSIZE*i+j ];
            }
        }
        if(compressible_char > 0 && compressible_char < DYNAMICWORDSIZE)
        {
            bool isSame = false;
            for(int j = 0; j < pos; j++)
            {
                if(extracted_patterns[j] == patterns_temp[i] || patterns_temp[i] == 51 || patterns_temp[i] == 3 || patterns_temp[i] == 240 || patterns_temp[i] == 15)
                {
                    isSame = true;
                    break;
                } 
            }
            if(!isSame)
            {
                compressible_chars[pos] = compressible_char;
                extracted_patterns[pos] = patterns_temp[i];
                DynamicPatterns[pos + FPCCOUNT + BDICOUNT] = pos + FPCCOUNT + BDICOUNT;
                CompressBytes[pos + FPCCOUNT + BDICOUNT] =   min_compression_counter * compressible_char / 2 - min_compression_counter * 3 / 8;//min_compression_counter * compressible_char / 2 - min_compression_counter * 3 / 8) / (DYNAMICWORDSIZE / 2);
                //min_compression_counter * (70 - 8 * compressible_char);
                //SampleCompressBytes[pos] = min_compression_counter * (70 - 8 * compressible_char);
                // std::cout<<"extracted_pattern: "<<extracted_patterns[pos]<<" CompressBytes["<<pos + FPCCOUNT + BDICOUNT<<"]: "<<CompressBytes[pos + FPCCOUNT + BDICOUNT]<<std::endl;
                pos++;
            }
            //(min_compression_counter * compressible_char / 2 - min_compression_counter * 3 / 8) / (DYNAMICWORDSIZE / 2);
        }
    }
    printf("\n");
    HeapSort(DynamicPatterns, CompressBytes, pos + FPCCOUNT + BDICOUNT, 4);
    int count;
    for(i = 0, count = pos + FPCCOUNT + BDICOUNT - 1; i < 4; i++, count--)
    {
        uint8_t pattern;
        uint64_t mask = 0;
        if(DynamicPatterns[count] < FPCCOUNT)
        {
            switch(DynamicPatterns[count])
            {
                case 0:
                    pattern = 3;
                    printf("3\t"); //00000011
                    for(j=0; j<8;j++)
                    {
                        if((pattern & 0x1) == 0)
                        {
                            mask = mask | (0xF << (j*4));
                        }
                        pattern = pattern >> 1;
                    }
                    masks[mask_pos] = mask;
                    // std::cout<<" mask: "<<mask<<" ";
                    //这时候Mask是0x ffffff00
                    compressibleChars[mask_pos++] = 6;
                    
                    break;
				case 1:
                    pattern = 51;
                    printf("51\t");//00110011
                    for(j=0; j<8;j++)
                    {
                        if((pattern & 0x1) == 0)
                        {
                            mask = mask | (0xF << (j*4));
                        }
                        pattern = pattern >> 1;
                    }
                    masks[mask_pos] = mask;
                    // std::cout<<" mask: "<<mask<<" ";
                    // 这时候mask是0xff00ff00
                    compressibleChars[mask_pos++] = 4;
                    
                    break;
                case 2:
                    special_pattern_flag[0] = true;
                    printf("sameBytes\t");
                    break;
            }
            
            //DynamicPatterns[count]
        }
        else if(DynamicPatterns[count] < FPCCOUNT + BDICOUNT)
        {
            special_pattern_flag[1+DynamicPatterns[count] - FPCCOUNT] = true;
            switch(DynamicPatterns[count] - FPCCOUNT)
            {
                case 0:
                    printf("same64bitsWord\t");
                    break;
                case 1:
                    printf("1-8\t");
                    break;
                case 2:
                    printf("2-8\t");
                    break;
                case 3:
                    printf("4-8\t");
                    break;
                case 4:
                    printf("same32bitsWord\t");
                    break;
                case 5:
                    printf("1-4\t");
                    break;
                case 6:
                    printf("2-4\t");
                    break;
                case 7:
                    printf("1-2\t");
                    break;
            }
            //DynamicPatterns[count] - FPCCOUNT
        }
        else
        {
			pattern = extracted_patterns[DynamicPatterns[count] - FPCCOUNT - BDICOUNT];
			printf("%d\t", pattern);
			for(j=0; j<8;j++)
			{
				if((pattern & 0x1) == 0)
				{
					mask = mask | (0xF << (j*4));
				}
				pattern = pattern >> 1;
			}
            // std::cout<<" mask: "<<mask<<" ";
            // std::cout<<" comChar: "<<compressible_chars[DynamicPatterns[count] - FPCCOUNT - BDICOUNT]<<" ";
            masks[mask_pos] = mask;
            compressibleChars[mask_pos++] = compressible_chars[DynamicPatterns[count] - FPCCOUNT - BDICOUNT];
            
            //extracted_patterns[DynamicPatterns[count] - FPCCOUNT - BDICOUNT]
        }
    }
    printf("\n");
    for(i = 0; i<mask_pos; i++)
    {
        // std::cout<<" mask: "<<masks[i]<<" comChar: "<<compressibleChars[i]<<std::endl;
    }
    sample_flag = false;
	return 1;
}

void FRFCFS::HeapAdjust(uint64_t pattern_array[], uint64_t bytes_array[],int pos,int nLength)
{
    int i, nChild;
    uint64_t nTemp;
    for(i = pos; 2*i+1 < nLength; i = nChild)
    {
        nChild = 2*i+1;
        if((nChild < nLength-1) && (bytes_array[nChild+1] > bytes_array[nChild]))
            ++nChild;
        if(bytes_array[nChild] > bytes_array[i])
        {
            nTemp = bytes_array[i];
            bytes_array[i] = bytes_array[nChild];
            bytes_array[nChild] = nTemp;
            
            nTemp = pattern_array[i];
            pattern_array[i] = pattern_array[nChild];
            pattern_array[nChild] = nTemp;
        }
        else break;
    }
}

void FRFCFS::HeapSort(uint64_t pattern_array[], uint64_t bytes_array[],int length, int topk)
{
    int i, count;
    for(i = length/2-1; i >= 0; --i)
        HeapAdjust(pattern_array, bytes_array, i , length);
    for(i = length-1, count = 0; count < topk; --i, count++)
    {
        pattern_array[i]=pattern_array[0]^pattern_array[i];
        pattern_array[0]=pattern_array[0]^pattern_array[i];
        pattern_array[i]=pattern_array[0]^pattern_array[i];
        
        bytes_array[i]=bytes_array[0]^bytes_array[i];
        bytes_array[0]=bytes_array[0]^bytes_array[i];
        bytes_array[i]=bytes_array[0]^bytes_array[i];
        HeapAdjust(pattern_array, bytes_array,0,i);
    }
}
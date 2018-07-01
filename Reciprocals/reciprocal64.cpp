#include <stdint.h>
#include <intrin.h>

static const uint32_t R[256] = {
	0xfe01be73u,0xfa0efe02u,0xfa118cdau,0xf249f625u,0xf630ce32u,0xeab1ee77u,0xf25f28d1u,0xe345e6f6u,
	0xee9c4562u,0xdc04dfa0u,0xeae7cedfu,0xd4edd874u,0xe7417330u,0xcdfed170u,0xe3a8e17fu,0xc736ca95u,
	0xe01dcd03u,0xc094c3e0u,0xdc9fea3bu,0xba16bd50u,0xd92eef95u,0xb3bdb6e5u,0xd5ca9688u,0xad86b09du,
	0xd2729981u,0xa771aa78u,0xcf26b64du,0xa17da473u,0xcbe6aae6u,0x9ba99e8fu,0xc8b237ffu,0x95f498cbu,
	0xc589204fu,0x905d9325u,0xc26b27abu,0x8ae38d9du,0xbf58143fu,0x85868831u,0xbc4fac82u,0x804582e2u,
	0xb951b97fu,0x7b1f7dafu,0xb65e05d1u,0x76137896u,0xb3745d33u,0x71217397u,0xb0948c50u,0x6c476eb1u,
	0xadbe61e8u,0x678669e4u,0xaaf1ada4u,0x62dd652fu,0xa82e40adu,0x5e4b6091u,0xa573ecafu,0x59cf5c0au,
	0xa2c2856du,0x55695799u,0xa019df1cu,0x5118533eu,0x9d79cfcdu,0x4cdd4ef8u,0x9ae22da2u,0x48b54ac6u,
	0x9852d02eu,0x44a146a9u,0x95cb90b5u,0x40a1429fu,0x934c4884u,0x3cb33ea8u,0x90d4d2a6u,0x38d83ac3u,
	0x8e65099fu,0x350f36f1u,0x8bfcca3du,0x31573331u,0x899bf212u,0x2db02f82u,0x87425f4du,0x2a1a2be3u,
	0x84eff042u,0x26942855u,0x82a48484u,0x231f24d8u,0x805ffd31u,0x1fb82169u,0x7e223a4du,0x1c611e0bu,
	0x7beb1f14u,0x19191abbu,0x79ba8d10u,0x15df177au,0x779067d1u,0x12b31447u,0x756c934eu,0xf951122u,
	0x734ef3ceu,0xc840e0bu,0x71376ed5u,0x9810b01u,0x6f25e9fbu,0x68b0804u,0x6d1a4b55u,0x3a10514u,
	0x6b1479e7u,0xc40231u,0x69145d70u,0xfdf2ff59u,0x6719dd5au,0xfb2dfc8eu,0x6524e2c2u,0xf873f9ceu,
	0x6335561bu,0xf5c4f71au,0x614b2163u,0xf320f471u,0x5f662e94u,0xf088f1d3u,0x5d86681fu,0xedfaef3fu,
	0x5babb887u,0xeb76ecb6u,0x59d60b2eu,0xe8fcea38u,0x58054c96u,0xe68ce7c3u,0x5639688fu,0xe426e558u,
	0x54724b96u,0xe1cae2f7u,0x52afe2feu,0xdf77e09fu,0x50f21bc9u,0xdd2dde51u,0x4f38e452u,0xdaecdc0bu,
	0x4d842a0eu,0xd8b3d9cfu,0x4bd3dc5au,0xd684d79bu,0x4a27e9dcu,0xd45dd56fu,0x4880415eu,0xd23ed34cu,
	0x46dcd2bau,0xd027d131u,0x453d8df4u,0xce18cf1eu,0x43a2630fu,0xcc11cd13u,0x420b4269u,0xca11cb10u,
	0x40781d76u,0xc819c914u,0x3ee8e49au,0xc628c720u,0x3d5d89d3u,0xc43fc533u,0x3bd5fe6du,0xc25cc34du,
	0x3a5234b6u,0xc081c16eu,0x38d21e75u,0xbeacbf96u,0x3755aec4u,0xbcdebdc4u,0x35dcd7e8u,0xbb16bbf9u,
	0x34678cbfu,0xb955ba35u,0x32f5c0f8u,0xb79ab877u,0x318767edu,0xb5e5b6bfu,0x301c7585u,0xb437b50du,
	0x2eb4dccau,0xb28eb362u,0x2d5092f9u,0xb0ebb1bcu,0x2bef8bfau,0xaf4eb01cu,0x2a91bc94u,0xadb7ae82u,
	0x2937195cu,0xac25acedu,0x27df970eu,0xaa98ab5eu,0x268b2ba0u,0xa911a9d4u,0x2539cbffu,0xa78fa84fu,
	0x23eb6d84u,0xa612a6d0u,0x22a006aau,0xa49aa555u,0x21578c79u,0xa327a3e0u,0x2011f59bu,0xa1b9a270u,
	0x1ecf388eu,0xa050a104u,0x1d8f4b94u,0x9eec9f9du,0x1c5224ebu,0x9d8c9e3bu,0x1b17bb87u,0x9c309cdeu,
	0x19e0073fu,0x9ada9b84u,0x18aafd9au,0x99879a30u,0x17789704u,0x983998e0u,0x1648caa2u,0x96ef9794u,
	0x151b9001u,0x95a9964cu,0x13f0deeau,0x94689508u,0x12c8aef3u,0x932a93c8u,0x11a2f800u,0x91f1928du,
	0x107fb258u,0x90bb9155u,0xf5ed663u,0x8f899021u,0xe405c34u,0x8e5b8ef1u,0xd243c18u,0x8d308dc5u,
	0xc0a6f45u,0x8c0a8c9du,0xaf2edf2u,0x8ae68b78u,0x9ddb19bu,0x89c78a56u,0x8cab299u,0x88ab8938u,
	0x7b9e9d5u,0x8792881eu,0x6ab5148u,0x867d8707u,0x59ee1bcu,0x856b85f3u,0x494945cu,0x845c84e3u,
	0x38c62ffu,0x835083d6u,0x286478bu,0x824882ccu,0x1823b26u,0x814281c5u,0x803807u,0x804080c1u,
};

struct SReciprocal
{
	uint32_t C0;
	uint16_t C1;
	uint16_t C2;
};

#define RCP ((SReciprocal*)R)

//uint64_t get_reciprocal(uint32_t a)
//{
//	unsigned long shift;
//	_BitScanReverse64(&shift, a);
//	shift = 31 - shift;
//	const uint64_t a1 = a << shift;
//
//	const uint32_t index1 = (static_cast<uint32_t>(a1 << 1) >> 25);
//	const int32_t index2 = static_cast<int32_t>(a1 & 16777215) - 8388607;
//	const int64_t r1 = RCP[index1].C0 | 0x100000000LL;
//	const int32_t r2_1 = RCP[index1].C1;
//	const int32_t r2_2 = RCP[index1].C2;
//	const int32_t r2 = ((index2 < 0) ? r2_1 : r2_2) | ((index1 < 53) ? 65536 : 0);
//	uint64_t r = (r1 - ((static_cast<int64_t>(r2) * index2) >> 15)) << shift;
//
//	uint64_t lo, hi, hi2;
//
//	lo = _umul128(r, a, &hi);
//	_subborrow_u64(_subborrow_u64(0, a >> 1, lo, &lo), 2, hi, &hi);
//	_umul128(r, lo, &hi2);
//	r = hi2 + (hi ? r : 0);
//
//	if (a <= 76944024)
//	{
//		lo = _umul128(r, a, &hi);
//		_subborrow_u64(_subborrow_u64(0, a >> 1, lo, &lo), 2, hi, &hi);
//		_umul128(r, lo, &hi2);
//		r = hi2 + (hi ? r : 0);
//	}
//
//	return r;
//}

// Simpler and faster version only for a >= 2^31
// This one makes fast_div 7% faster than DIV instruction
uint64_t get_reciprocal(uint32_t a)
{
	const uint32_t index1 = (static_cast<uint32_t>(a << 1) >> 25);
	const int32_t index2 = static_cast<int32_t>(a & 16777215) - 8388607;
	const int64_t r1 = RCP[index1].C0 | 0x100000000LL;
	const int32_t r2_1 = RCP[index1].C1;
	const int32_t r2_2 = RCP[index1].C2;
	const int32_t r2 = ((index2 < 0) ? r2_1 : r2_2) | ((index1 < 53) ? 65536 : 0);
	uint64_t r = r1 - ((static_cast<int64_t>(r2) * index2) >> 15);

	uint64_t lo, hi, hi2;

	lo = _umul128(r, a, &hi);
	_subborrow_u64(_subborrow_u64(0, a >> 1, lo, &lo), 2, hi, &hi);
	_umul128(r, lo, &hi2);
	r = hi2 + (hi ? r : 0);

	return r;
}

void fast_div(uint64_t a, uint32_t b, uint64_t &q, uint64_t &r)
{
	_umul128(a, get_reciprocal(b), &q);
	int64_t tmp = a - q * b;
	if (tmp < 0)
	{
		--q;
		tmp += b;
	}
	else if (tmp >= b)
	{
		++q;
		tmp -= b;
	}
	r = static_cast<uint64_t>(tmp);
}
